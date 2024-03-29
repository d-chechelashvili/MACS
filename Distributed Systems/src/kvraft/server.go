package kvraft

import (
	"log"
	"sync"
	"sync/atomic"
	"time"

	"../labgob"
	"../labrpc"
	"../raft"
)

const Debug = 0

func DPrintf(format string, a ...interface{}) (n int, err error) {
	if Debug > 0 {
		log.Printf(format, a...)
	}
	return
}

type Op struct {
	// Your definitions here.
	// Field names must start with capital letters,
	// otherwise RPC will break.
	Type      string
	Key       string
	Value     string
	ClerkId   int64
	CommandId int
}

type ProcessedRequest struct {
	clerkId   int64
	commandId int
	answer    string
}

type KVServer struct {
	mu      sync.Mutex
	me      int
	rf      *raft.Raft
	applyCh chan raft.ApplyMsg
	dead    int32 // set by Kill()

	maxraftstate int // snapshot if log grows this big

	// Your definitions here.
	kvStore              map[string]string
	lastProcessedRequest map[int64]ProcessedRequest
	appliedRequests      map[int]ProcessedRequest
}

func (kv *KVServer) Get(args *GetArgs, reply *GetReply) {
	// Your code here.
	_, isLeader := kv.rf.GetState()
	if !isLeader {
		reply.Err = ErrWrongLeader
		return
	}

	op := Op{"Get", args.Key, "", args.ClerkId, args.CommandId}
	index, _, _ := kv.rf.Start(op)
	time.Sleep(520 * time.Millisecond)

	kv.mu.Lock()
	defer kv.mu.Unlock()
	if req, found := kv.appliedRequests[index]; found {
		if req.clerkId == op.ClerkId && req.commandId == op.CommandId {
			reply.Err = OK
			reply.Value = req.answer
			return
		}
	}

	reply.Err = ErrNoKey
}

func (kv *KVServer) PutAppend(args *PutAppendArgs, reply *PutAppendReply) {
	// Your code here.
	_, isLeader := kv.rf.GetState()

	if !isLeader {
		reply.Err = ErrWrongLeader
		return
	}

	op := Op{args.Op, args.Key, args.Value, args.ClerkId, args.CommandId}
	index, _, _ := kv.rf.Start(op)
	time.Sleep(520 * time.Millisecond)

	kv.mu.Lock()
	defer kv.mu.Unlock()
	if req, found := kv.appliedRequests[index]; found {
		if req.clerkId == args.ClerkId && req.commandId == args.CommandId {
			reply.Err = OK
			return
		}
	}

	reply.Err = ErrNoKey
}

//
// the tester calls Kill() when a KVServer instance won't
// be needed again. for your convenience, we supply
// code to set rf.dead (without needing a lock),
// and a killed() method to test rf.dead in
// long-running loops. you can also add your own
// code to Kill(). you're not required to do anything
// about this, but it may be convenient (for example)
// to suppress debug output from a Kill()ed instance.
//
func (kv *KVServer) Kill() {
	atomic.StoreInt32(&kv.dead, 1)
	kv.rf.Kill()
	// Your code here, if desired.
}

func (kv *KVServer) killed() bool {
	z := atomic.LoadInt32(&kv.dead)
	return z == 1
}

func (kv *KVServer) checkExecuted(clerkId int64, cmdId int) (bool, string) {
	kv.mu.Lock()
	defer kv.mu.Unlock()
	if ans, found := kv.lastProcessedRequest[clerkId]; found {
		if ans.commandId == cmdId {
			return true, ans.answer
		}
	}
	return false, ""
}

func (kv *KVServer) applyListener() {
	for !kv.killed() {
		applyMsg := <-kv.applyCh
		op := applyMsg.Command.(Op)

		if executed, _ := kv.checkExecuted(op.ClerkId, op.CommandId); !executed {
			if op.Type == "Put" {
				kv.kvStore[op.Key] = op.Value
			} else if op.Type == "Append" {
				kv.kvStore[op.Key] += op.Value
			}
		}

		kv.mu.Lock()
		processedRequest := ProcessedRequest{clerkId: op.ClerkId, commandId: op.CommandId, answer: kv.kvStore[op.Key]}
		kv.appliedRequests[applyMsg.CommandIndex] = processedRequest
		if op.Type != "Get" {
			kv.lastProcessedRequest[op.ClerkId] = processedRequest
		}
		kv.mu.Unlock()
	}
}

//
// servers[] contains the ports of the set of
// servers that will cooperate via Raft to
// form the fault-tolerant key/value service.
// me is the index of the current server in servers[].
// the k/v server should store snapshots through the underlying Raft
// implementation, which should call persister.SaveStateAndSnapshot() to
// atomically save the Raft state along with the snapshot.
// the k/v server should snapshot when Raft's saved state exceeds maxraftstate bytes,
// in order to allow Raft to garbage-collect its log. if maxraftstate is -1,
// you don't need to snapshot.
// StartKVServer() must return quickly, so it should start goroutines
// for any long-running work.
//
func StartKVServer(servers []*labrpc.ClientEnd, me int, persister *raft.Persister, maxraftstate int) *KVServer {
	// call labgob.Register on structures you want
	// Go's RPC library to marshall/unmarshall.
	labgob.Register(Op{})

	kv := new(KVServer)
	kv.me = me
	kv.maxraftstate = maxraftstate

	// You may need initialization code here.
	kv.applyCh = make(chan raft.ApplyMsg)
	kv.rf = raft.Make(servers, me, persister, kv.applyCh)

	kv.kvStore = make(map[string]string)
	kv.appliedRequests = make(map[int]ProcessedRequest)
	kv.lastProcessedRequest = make(map[int64]ProcessedRequest)

	// You may need initialization code here.
	go kv.applyListener()

	return kv
}
