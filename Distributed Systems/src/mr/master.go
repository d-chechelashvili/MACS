package mr

import (
	"log"
	"strconv"
	"sync"
	"time"
)
import "net"
import "os"
import "net/rpc"
import "net/http"

const (
	UNSTARTED = iota
	RUNNING
	COMPLETED
)

const (
	MAP_TYPE = iota
	REDUCE_TYPE
	WAIT_TYPE
	DONE_TYPE
)

type TaskInfo struct {
	Status    int
	StartTime time.Time
}

type Master struct {
	// Your definitions here.
	mapTasks    map[string]TaskInfo // filename, status
	nReduce     int                 // number of reduce tasks
	reduceTasks map[string]TaskInfo // task number, status
	mapTaskNum  map[string]int
	lock        sync.Mutex
}

// Your code here -- RPC handlers for the worker to call.
func (m *Master) GetTask(args *TaskCallArgs, reply *TaskCallReply) error {
	m.lock.Lock()
	defer m.lock.Unlock()

	// map phase
	mapTasksDone, mapTaskIssued := checkMapTasks(m, reply)
	if mapTaskIssued {return nil}
	if !mapTasksDone {
		// if task was not issued but also map tasks aren't done worker should wait
		reply.TaskType = WAIT_TYPE
		return nil
	}

	// reduce phase
	reduceTasksDone, reduceTaskIssued := checkReduceTasks(m, reply)
	if reduceTaskIssued {return nil}
	if !reduceTasksDone {
		// if task was not issued but also reduce tasks aren't done worker should wait
		reply.TaskType = WAIT_TYPE
		return nil
	}

	// if no task was issued and all tasks are done then worker is done too
	reply.TaskType = DONE_TYPE
	return nil
}

func checkReduceTasks(m *Master, reply *TaskCallReply) (bool, bool) {
	doneReduceTasks := 0
	for taskNumber, info := range m.reduceTasks {
		if info.Status == RUNNING && time.Since(info.StartTime) >= 10*time.Second {
			m.reduceTasks[taskNumber] = TaskInfo{UNSTARTED, time.Now()}
		}

		if info.Status == UNSTARTED {
			m.reduceTasks[taskNumber] = TaskInfo{RUNNING, time.Now()}
			reply.TaskName = taskNumber
			reply.TaskNumber, _ = strconv.Atoi(taskNumber)
			reply.TaskType = REDUCE_TYPE
			reply.N = len(m.mapTaskNum)
			return false, true
		} else if info.Status == COMPLETED {
			doneReduceTasks++
		}
	}

	return doneReduceTasks == len(m.reduceTasks), false
}

func checkMapTasks(m *Master, reply *TaskCallReply) (bool, bool) {
	doneMapTasks := 0
	for filename, info := range m.mapTasks {
		if info.Status == RUNNING && time.Since(info.StartTime) >= 10*time.Second {
			m.mapTasks[filename] = TaskInfo{UNSTARTED, time.Now()}
		}

		if info.Status == UNSTARTED {
			m.mapTasks[filename] = TaskInfo{RUNNING, time.Now()}
			reply.TaskName = filename
			reply.TaskNumber = m.mapTaskNum[filename]
			reply.TaskType = MAP_TYPE
			reply.N = m.nReduce
			return false, true
		} else if info.Status == COMPLETED {
			doneMapTasks++
		}
	}

	return doneMapTasks == len(m.mapTasks), false
}

func (m *Master) TaskDone(args *TaskCallReply, reply *TaskCallReply) error {
	m.lock.Lock()
	defer m.lock.Unlock()
	if args.TaskType == MAP_TYPE {
		info := m.mapTasks[args.TaskName]
		info.Status = COMPLETED
		m.mapTasks[args.TaskName] = info
	} else if args.TaskType == REDUCE_TYPE {
		info := m.reduceTasks[args.TaskName]
		info.Status = COMPLETED
		m.reduceTasks[args.TaskName] = info
	}
	return nil
}

//
// an example RPC handler.
//
// the RPC argument and reply types are defined in rpc.go.
//
func (m *Master) Example(args *ExampleArgs, reply *ExampleReply) error {
	reply.Y = args.X + 1
	return nil
}

//
// start a thread that listens for RPCs from worker.go
//
func (m *Master) server() {
	rpc.Register(m)
	rpc.HandleHTTP()
	//l, e := net.Listen("tcp", ":1234")
	sockname := masterSock()
	os.Remove(sockname)
	l, e := net.Listen("unix", sockname)
	if e != nil {
		log.Fatal("listen error:", e)
	}
	go http.Serve(l, nil)
}

//
// main/mrmaster.go calls Done() periodically to find out
// if the entire job has finished.
//
func (m *Master) Done() bool {
	//ret := false
	// Your code here.
	for _, info := range m.reduceTasks {
		if info.Status != COMPLETED {
			return false
		}
	}
	return true
}

//
// create a Master.
// main/mrmaster.go calls this function.
// nReduce is the number of reduce tasks to use.
//
func MakeMaster(files []string, nReduce int) *Master {
	m := Master{}
	// Your code here.
	m.mapTasks = make(map[string]TaskInfo)
	m.reduceTasks = make(map[string]TaskInfo)
	m.nReduce = nReduce
	m.mapTaskNum = make(map[string]int)

	for i, filename := range files {
		m.mapTasks[filename] = TaskInfo{UNSTARTED, time.Time{}}
		m.mapTaskNum[filename] = i
	}

	for i := 0; i < nReduce; i++ {
		m.reduceTasks[strconv.Itoa(i)] = TaskInfo{UNSTARTED, time.Time{}}
	}
	m.server()
	return &m
}
