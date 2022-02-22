package raft

//
// this is an outline of the API that raft must expose to
// the service (or tester). see comments below for
// each of these functions for more details.
//
// rf = Make(...)
//   create a new Raft server.
// rf.Start(command interface{}) (index, term, isleader)
//   start agreement on a new log entry
// rf.GetState() (term, isLeader)
//   ask a Raft for its current term, and whether it thinks it is leader
// ApplyMsg
//   each time a new entry is committed to the log, each Raft peer
//   should send an ApplyMsg to the service (or tester)
//   in the same server.
//

import (
	"bytes"
	"math/rand"
	"sync"
	"sync/atomic"
	"time"

	"../labgob"
	"../labrpc"
)

const (
	FOLLOWER = iota
	CANDIDATE
	LEADER
)

const (
	HEARTBEAT_TIMEOUT_RAND_RANGE              = 350
	HEARTBEAT_MIN_TIMEOUT                     = 400
	CHECKING_SLEEP_INTERVAL                   = 80
	LEADER_ADDITIONAL_CHECKING_SLEEP_INTERVAL = 35
	APPLIER_SLEEP_INTERVAL                    = 13
)

//
// as each Raft peer becomes aware that successive log entries are
// committed, the peer should send an ApplyMsg to the service (or
// tester) on the same server, via the applyCh passed to Make(). set
// CommandValid to true to indicate that the ApplyMsg contains a newly
// committed log entry.
//
// in Lab 3 you'll want to send other kinds of messages (e.g.,
// snapshots) on the applyCh; at that point you can add fields to
// ApplyMsg, but set CommandValid to false for these other uses.
//
type ApplyMsg struct {
	CommandValid bool
	Command      interface{}
	CommandIndex int
}

type Log struct {
	Command interface{}
	Term    int
}

//
// A Go object implementing a single Raft peer.
//

type Raft struct {
	mu        sync.Mutex          // Lock to protect shared access to this peer's state
	peers     []*labrpc.ClientEnd // RPC end points of all peers
	persister *Persister          // Object to hold this peer's persisted state
	me        int                 // this peer's index into peers[]
	dead      int32               // set by Kill()

	// Your data here (2A, 2B, 2C).
	// Look at the paper's Figure 2 for a description of what
	// state a Raft server must maintain.

	//my Data
	status            int
	lastHeartbeatTime time.Time
	voteCount         int
	applyCh           chan ApplyMsg

	//data from Figure 2
	//Persistent state on all servers:
	currentTerm int   //latest term server has seen
	votedFor    int   //candidateId that received vote in current term (or -1 if none)
	log         []Log //log entries; each entry contains command for state machine, and term when
	//entry was received by leader (first index is 1)

	//Volatile state on all servers:
	commitIndex int //index of highest log entry known to be committed
	lastApplied int //index of highest log entry applied to state machine

	//Volatile state on leaders:
	nextIndex  []int //for each server, index of the next log entry to send to that server
	matchIndex []int //for each server, index of highest log entry known to be replicated on server
}

// return currentTerm and whether this server
// believes it is the leader.
func (rf *Raft) GetState() (int, bool) {
	var term int
	var isleader bool
	// Your code here (2A).
	rf.mu.Lock()
	defer rf.mu.Unlock()
	term = rf.currentTerm
	isleader = (rf.status == LEADER)
	return term, isleader
}

//
// save Raft's persistent state to stable storage,
// where it can later be retrieved after a crash and restart.
// see paper's Figure 2 for a description of what should be persistent.
//
func (rf *Raft) persist() {
	// Your code here (2C).
	w := new(bytes.Buffer)
	e := labgob.NewEncoder(w)
	e.Encode(rf.currentTerm)
	e.Encode(rf.votedFor)
	e.Encode(rf.log)
	data := w.Bytes()
	rf.persister.SaveRaftState(data)
}

//
// restore previously persisted state.
//
func (rf *Raft) readPersist(data []byte) {
	if data == nil || len(data) < 1 { // bootstrap without any state?
		return
	}
	// Your code here (2C).
	r := bytes.NewBuffer(data)
	d := labgob.NewDecoder(r)
	var currentTerm int
	var votedFor int
	var logs []Log
	if d.Decode(&currentTerm) != nil || d.Decode(&votedFor) != nil || d.Decode(&logs) != nil {
	} else {
		rf.currentTerm = currentTerm
		rf.votedFor = votedFor
		rf.log = logs
	}
}

//
// example RequestVote RPC arguments structure.
// field names must start with capital letters!
//
type RequestVoteArgs struct {
	// Your data here (2A, 2B).
	Term         int
	CandidateId  int
	LastLogIndex int
	LastLogTerm  int
}

//
// example RequestVote RPC reply structure.
// field names must start with capital letters!
//
type RequestVoteReply struct {
	// Your data here (2A).
	Term        int
	VoteGranted bool
}

//
// example RequestVote RPC handler.
//
func (rf *Raft) RequestVote(args *RequestVoteArgs, reply *RequestVoteReply) {
	// Your code here (2A, 2B).
	rf.mu.Lock()
	defer rf.mu.Unlock()

	reply.Term = rf.currentTerm

	if args.Term < rf.currentTerm {
		reply.VoteGranted = false
		return
	}

	if args.Term > rf.currentTerm {
		rf.convertToFollower(args.Term, -1)
	}

	//if haven't voted in this term and also requesters log is up to date to mine I can vote
	if (rf.votedFor < 0 || rf.votedFor == args.CandidateId) && logIsUpToDate(args, rf) {
		reply.VoteGranted = true
		rf.votedFor = args.CandidateId
		rf.persist()
		rf.lastHeartbeatTime = time.Now()
	} else {
		reply.VoteGranted = false
	}
}

// first compare terms, then compare indices
func logIsUpToDate(args *RequestVoteArgs, rf *Raft) bool {
	lastLog := rf.log[len(rf.log)-1]
	if args.LastLogTerm < lastLog.Term {
		return false
	}

	if args.LastLogTerm > lastLog.Term {
		return true
	}
	return args.LastLogIndex >= (len(rf.log) - 1)
}

func (rf *Raft) convertToFollower(term int, votedFor int) {
	defer rf.persist()
	rf.currentTerm = term
	rf.status = FOLLOWER
	rf.votedFor = votedFor
	rf.voteCount = 0
}

//
// example code to send a RequestVote RPC to a server.
// server is the index of the target server in rf.peers[].
// expects RPC arguments in args.
// fills in *reply with RPC reply, so caller should
// pass &reply.
// the types of the args and reply passed to Call() must be
// the same as the types of the arguments declared in the
// handler function (including whether they are pointers).
//
// The labrpc package simulates a lossy network, in which servers
// may be unreachable, and in which requests and replies may be lost.
// Call() sends a request and waits for a reply. If a reply arrives
// within a timeout interval, Call() returns true; otherwise
// Call() returns false. Thus Call() may not return for a while.
// A false return can be caused by a dead server, a live server that
// can't be reached, a lost request, or a lost reply.
//
// Call() is guaranteed to return (perhaps after a delay) *except* if the
// handler function on the server side does not return.  Thus there
// is no need to implement your own timeouts around Call().
//
// look at the comments in ../labrpc/labrpc.go for more details.
//
// if you're having trouble getting RPC to work, check that you've
// capitalized all field names in structs passed over RPC, and
// that the caller passes the address of the reply struct with &, not
// the struct itself.
//
func (rf *Raft) sendRequestVote(server int, args *RequestVoteArgs, reply *RequestVoteReply) bool {
	ok := rf.peers[server].Call("Raft.RequestVote", args, reply)
	return ok
}

//start threads for each peer to ask for votes
func (rf *Raft) sendRequestVotes(args *RequestVoteArgs) {
	for peer := range rf.peers {
		if peer != rf.me {
			go rf.requestVoteSender(peer, args)
		}
	}
}

//asks for vote to designated peer and takes action depending on reply
func (rf *Raft) requestVoteSender(peer int, args *RequestVoteArgs) {
	reply := RequestVoteReply{}

	rf.mu.Lock()
	if rf.status != CANDIDATE {
		rf.mu.Unlock()
		return
	}
	rf.mu.Unlock()

	rf.sendRequestVote(peer, args, &reply)

	rf.mu.Lock()
	defer rf.mu.Unlock()

	if reply.Term > rf.currentTerm {
		rf.convertToFollower(reply.Term, -1)
		return
	}

	if rf.currentTerm != args.Term { //If term changed before answer was returned discard answer
		return
	}

	if reply.VoteGranted {
		rf.voteCount++
	}
}

type AppendEntriesArgs struct {
	Term         int
	LeaderId     int
	PrevLogIndex int
	PrevLogTerm  int
	Entries      []Log
	LeaderCommit int
}

type AppendEntriesReply struct {
	Term          int
	Success       bool
	ConflictIndex int
	ConflictTerm  int
}

// AppendEntries RPC handler
func (rf *Raft) AppendEntries(args *AppendEntriesArgs, reply *AppendEntriesReply) {
	rf.mu.Lock()
	defer rf.mu.Unlock()

	reply.Term = rf.currentTerm

	if args.Term < rf.currentTerm {
		reply.Success = false
		return
	}

	//if term is >= then this AppendEntry is from leader
	if args.Term > rf.currentTerm {
		rf.convertToFollower(args.Term, -1)
	} else {
		rf.convertToFollower(args.Term, rf.votedFor)
	}
	rf.lastHeartbeatTime = time.Now()

	//If we don't even have PrevLog in our Log return false
	if args.PrevLogIndex >= len(rf.log) {
		reply.Success = false
		reply.ConflictIndex = len(rf.log)
		reply.ConflictTerm = -2
		return
	}

	//If there is conflict with prevLogs return false and truncate our log
	if rf.log[args.PrevLogIndex].Term != args.PrevLogTerm {
		reply.Success = false
		reply.ConflictTerm = rf.log[args.PrevLogIndex].Term
		//find conflict index
		for reply.ConflictIndex = args.PrevLogIndex; reply.ConflictIndex > 0; {
			if rf.log[reply.ConflictIndex-1].Term != reply.ConflictTerm {
				break
			}
			reply.ConflictIndex--
		}
		rf.log = rf.log[:args.PrevLogIndex]
		defer rf.persist()
		return
	}

	rf.updateLog(args)
	rf.persist()
	reply.Success = true

	//updates our commitIndex
	if args.LeaderCommit > rf.commitIndex {
		indexOfLastNewEntry := args.PrevLogIndex + len(args.Entries)
		rf.commitIndex = max(rf.commitIndex, min(args.LeaderCommit, indexOfLastNewEntry))
	}
}

//updates Log based on args.Entries
func (rf *Raft) updateLog(args *AppendEntriesArgs) {
	for i := 0; i < len(args.Entries); i++ {
		checkIndex := args.PrevLogIndex + 1 + i
		//if there is conflict we truncate and add new entries
		if checkIndex >= len(rf.log) || rf.log[checkIndex].Term != args.Entries[i].Term {
			rf.log = append(rf.log[:checkIndex], args.Entries[i:]...)
			break
		} else { // if no conflict we update only given entries and avoid truncating
			rf.log[checkIndex] = args.Entries[i]
		}
	}
}

func (rf *Raft) sendAppendEntries(server int, args *AppendEntriesArgs, reply *AppendEntriesReply) bool {
	ok := rf.peers[server].Call("Raft.AppendEntries", args, reply)
	return ok
}

//start threads for each peer to send AppendEntries
func (rf *Raft) sendHeartbeats() {
	for peer := range rf.peers {
		if peer != rf.me {
			args := AppendEntriesArgs{Term: rf.currentTerm, LeaderId: rf.me}
			args.Entries = append(args.Entries, rf.log[rf.nextIndex[peer]:]...)
			args.PrevLogIndex = rf.nextIndex[peer] - 1
			args.PrevLogTerm = rf.log[args.PrevLogIndex].Term
			args.LeaderCommit = rf.commitIndex
			go rf.heartbeatSender(peer, &args)
		}
	}
}

//sends AppendEntries to designated peer and takes action depending on reply
func (rf *Raft) heartbeatSender(peer int, args *AppendEntriesArgs) {
	reply := AppendEntriesReply{}

	rf.mu.Lock()
	if rf.status != LEADER {
		rf.mu.Unlock()
		return
	}
	rf.mu.Unlock()

	rf.sendAppendEntries(peer, args, &reply)

	rf.mu.Lock()
	defer rf.mu.Unlock()

	if reply.Term > rf.currentTerm {
		rf.convertToFollower(reply.Term, -1)
		return
	}
	if rf.currentTerm != args.Term { //If term changed before answer was returned discard answer
		return
	}

	if reply.Success {
		newMatchIndex := args.PrevLogIndex + len(args.Entries)
		rf.matchIndex[peer] = max(rf.matchIndex[peer], newMatchIndex)
		rf.nextIndex[peer] = rf.matchIndex[peer] + 1
	} else if reply.ConflictIndex > 0 { //If AppendEntries fails because of log inconsistency
		for i := len(rf.log) - 1; i > 0; i-- {
			if rf.log[i].Term == reply.ConflictTerm {
				rf.nextIndex[peer] = max(rf.matchIndex[peer]+1, min(rf.nextIndex[peer], i+1)) //we might receive outdated packets even from this term
				return
			}
		}
		rf.nextIndex[peer] = max(rf.matchIndex[peer]+1, min(rf.nextIndex[peer], reply.ConflictIndex)) //we might receive outdated packet even from this term
	}
}

//
// the service using Raft (e.g. a k/v server) wants to start
// agreement on the next command to be appended to Raft's log. if this
// server isn't the leader, returns false. otherwise start the
// agreement and return immediately. there is no guarantee that this
// command will ever be committed to the Raft log, since the leader
// may fail or lose an election. even if the Raft instance has been killed,
// this function should return gracefully.
//
// the first return value is the index that the command will appear at
// if it's ever committed. the second return value is the current
// term. the third return value is true if this server believes it is
// the leader.
//
func (rf *Raft) Start(command interface{}) (int, int, bool) {
	index := -1
	term := -1
	isLeader := true

	// Your code here (2B).
	rf.mu.Lock()
	defer rf.mu.Unlock()
	term = rf.currentTerm
	isLeader = (rf.status == LEADER)

	if isLeader {
		// add command to log
		index = len(rf.log)
		rf.log = append(rf.log, Log{Command: command, Term: term})
		rf.matchIndex[rf.me] = index
		rf.nextIndex[rf.me] = rf.matchIndex[rf.me] + 1
		defer rf.persist()
	}

	return index, term, isLeader
}

//
// the tester doesn't halt goroutines created by Raft after each test,
// but it does call the Kill() method. your code can use killed() to
// check whether Kill() has been called. the use of atomic avoids the
// need for a lock.
//
// the issue is that long-running goroutines use memory and may chew
// up CPU time, perhaps causing later tests to fail and generating
// confusing debug output. any goroutine with a long-running loop
// should call killed() to check whether it should stop.
//
func (rf *Raft) Kill() {
	atomic.StoreInt32(&rf.dead, 1)
	// Your code here, if desired.
}

func (rf *Raft) killed() bool {
	z := atomic.LoadInt32(&rf.dead)
	return z == 1
}

//leader routine work
func (rf *Raft) checkAsLeader() {
	rf.mu.Lock()
	defer rf.mu.Unlock()
	if rf.status != LEADER {
		return
	}
	rf.sendHeartbeats()
}

//checks if commitIndex can be updated and updates it if possible
func (rf *Raft) updateCommits() {
	for n := len(rf.log) - 1; n > rf.commitIndex; n-- {
		count := 0
		if rf.log[n].Term == rf.currentTerm { //to not commit logs from previous terms
			for peer := range rf.peers {
				if rf.matchIndex[peer] >= n {
					count++
				}
			}
		}

		if count > len(rf.peers)/2 { //if more than n/2 servers had this log replicated we can commit
			rf.commitIndex = n
			break
		}
	}
}

//candidate routine work
func (rf *Raft) checkAsCandidate() {
	rf.mu.Lock()
	defer rf.mu.Unlock()
	if rf.status != CANDIDATE {
		return
	}

	if rf.voteCount > len(rf.peers)/2 {
		//convert to Leader
		rf.status = LEADER
		rf.voteCount = 1
		rf.votedFor = rf.me
		rf.persist()
		//reinitialize nextIndex, matchIndex
		for i := 0; i < len(rf.peers); i++ {
			rf.nextIndex[i] = len(rf.log)
			rf.matchIndex[i] = 0
		}
		rf.sendHeartbeats() //send initial heartbeats
		return
	}

	electionTimeout := HEARTBEAT_MIN_TIMEOUT + rand.Intn(HEARTBEAT_TIMEOUT_RAND_RANGE)
	if time.Since(rf.lastHeartbeatTime) >= time.Duration(electionTimeout)*time.Millisecond {
		rf.startElection()
	}
}

func (rf *Raft) startElection() {
	rf.lastHeartbeatTime = time.Now()
	rf.currentTerm++
	rf.voteCount = 1
	rf.votedFor = rf.me
	rf.persist()
	args := RequestVoteArgs{CandidateId: rf.me, Term: rf.currentTerm,
		LastLogIndex: len(rf.log) - 1, LastLogTerm: rf.log[len(rf.log)-1].Term}

	rf.sendRequestVotes(&args)
}

//follower routine work
func (rf *Raft) checkAsFollower() {
	rf.mu.Lock()
	defer rf.mu.Unlock()

	// status might have been modified in the checking process
	if rf.status != FOLLOWER {
		return
	}

	electionTimeout := HEARTBEAT_MIN_TIMEOUT + rand.Intn(HEARTBEAT_TIMEOUT_RAND_RANGE)
	if time.Since(rf.lastHeartbeatTime) >= time.Duration(electionTimeout)*time.Millisecond {
		rf.status = CANDIDATE
		rf.startElection()
	}
}

// do routine work according to status and sleep
func (rf *Raft) checkTasks() {
	for !rf.killed() {
		rf.mu.Lock()
		status := rf.status
		rf.mu.Unlock()

		if status == FOLLOWER {
			rf.checkAsFollower()
		}

		if status == CANDIDATE {
			rf.checkAsCandidate()
		}

		if status == LEADER {
			rf.checkAsLeader()
			time.Sleep(LEADER_ADDITIONAL_CHECKING_SLEEP_INTERVAL * time.Millisecond)
		}

		time.Sleep(CHECKING_SLEEP_INTERVAL * time.Millisecond)
	}
}

//periodically checks if logs can be commited or applied and takes action accordingly
func (rf *Raft) applier() {
	for !rf.killed() {
		rf.mu.Lock()
		if rf.status == LEADER {
			rf.updateCommits()
		}
		for rf.lastApplied < rf.commitIndex {
			rf.lastApplied++
			rf.applyCh <- ApplyMsg{true, rf.log[rf.lastApplied].Command, rf.lastApplied}
		}
		rf.mu.Unlock()
		time.Sleep(APPLIER_SLEEP_INTERVAL * time.Millisecond)
	}
}

//
// the service or tester wants to create a Raft server. the ports
// of all the Raft servers (including this one) are in peers[]. this
// server's port is peers[me]. all the servers' peers[] arrays
// have the same order. persister is a place for this server to
// save its persistent state, and also initially holds the most
// recent saved state, if any. applyCh is a channel on which the
// tester or service expects Raft to send ApplyMsg messages.
// Make() must return quickly, so it should start goroutines
// for any long-running work.
//
func Make(peers []*labrpc.ClientEnd, me int,
	persister *Persister, applyCh chan ApplyMsg) *Raft {
	rf := &Raft{}
	rf.peers = peers
	rf.persister = persister
	rf.me = me

	// Your initialization code here (2A, 2B, 2C).
	rf.currentTerm = 0
	rf.lastApplied = 0
	rf.commitIndex = 0
	rf.votedFor = -1
	rf.voteCount = 0
	rf.log = make([]Log, 1)
	rf.log[0].Term = -1
	rf.lastHeartbeatTime = time.Now()
	rf.status = FOLLOWER
	rf.matchIndex = make([]int, len(peers))
	rf.nextIndex = make([]int, len(peers))
	rf.applyCh = applyCh

	// initialize from state persisted before a crash
	rf.readPersist(persister.ReadRaftState())

	go rf.checkTasks()
	go rf.applier()

	return rf
}

func min(a int, b int) int {
	if b < a {
		return b
	}
	return a
}

func max(a int, b int) int {
	if b > a {
		return b
	}
	return a
}
