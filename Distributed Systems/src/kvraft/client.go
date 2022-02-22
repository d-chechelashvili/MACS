package kvraft

import (
	"crypto/rand"
	"math/big"
	"sync"

	"../labrpc"
)

type Clerk struct {
	servers []*labrpc.ClientEnd
	// You will have to modify this struct.
	commandId int
	leader    int
	id        int64
	mu        sync.Mutex
}

func nrand() int64 {
	max := big.NewInt(int64(1) << 62)
	bigx, _ := rand.Int(rand.Reader, max)
	x := bigx.Int64()
	return x
}

func MakeClerk(servers []*labrpc.ClientEnd) *Clerk {
	ck := new(Clerk)
	ck.servers = servers
	// You'll have to add code here.
	ck.commandId = 0
	ck.leader = int(nrand() % int64(len(servers)))
	ck.id = nrand()
	return ck
}

//
// fetch the current value for a key.
// returns "" if the key does not exist.
// keeps trying forever in the face of all other errors.
//
// you can send an RPC with code like this:
// ok := ck.servers[i].Call("KVServer.Get", &args, &reply)
//
// the types of args and reply (including whether they are pointers)
// must match the declared types of the RPC handler function's
// arguments. and reply must be passed as a pointer.
//
func (ck *Clerk) Get(key string) string {
	// You will have to modify this function.
	ck.mu.Lock()
	leader := ck.leader
	ck.commandId++
	args := GetArgs{Key: key, CommandId: ck.commandId, ClerkId: ck.id}
	ck.mu.Unlock()

	value := ""
	for {
		reply := GetReply{}
		ok := ck.servers[leader].Call("KVServer.Get", &args, &reply)

		if !ok || reply.Err == ErrWrongLeader {
			leader = (leader + 1) % len(ck.servers)
		} else {
			ck.mu.Lock()
			ck.leader = leader
			ck.mu.Unlock()
			if reply.Err == OK {
				value = reply.Value
				break
			}
		}
	}
	return value
}

//
// shared by Put and Append.
//
// you can send an RPC with code like this:
// ok := ck.servers[i].Call("KVServer.PutAppend", &args, &reply)
//
// the types of args and reply (including whether they are pointers)
// must match the declared types of the RPC handler function's
// arguments. and reply must be passed as a pointer.
//
func (ck *Clerk) PutAppend(key string, value string, op string) {
	// You will have to modify this function.
	ck.mu.Lock()
	leader := ck.leader
	ck.commandId++
	args := PutAppendArgs{Key: key, Value: value, Op: op, CommandId: ck.commandId, ClerkId: ck.id}
	ck.mu.Unlock()

	for {
		reply := PutAppendReply{}
		ok := ck.servers[leader].Call("KVServer.PutAppend", &args, &reply)

		if !ok || reply.Err == ErrWrongLeader {
			leader = (leader + 1) % len(ck.servers)
		} else {
			ck.mu.Lock()
			ck.leader = leader
			ck.mu.Unlock()
			if reply.Err == OK {
				break
			}
		}
	}
}

func (ck *Clerk) Put(key string, value string) {
	ck.PutAppend(key, value, "Put")
}
func (ck *Clerk) Append(key string, value string) {
	ck.PutAppend(key, value, "Append")
}
