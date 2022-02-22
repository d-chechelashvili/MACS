package mr

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"os"
	"strconv"
	"time"
)
import "log"
import "net/rpc"
import "hash/fnv"

//
// Map functions return a slice of KeyValue.
//
type KeyValue struct {
	Key   string
	Value string
}

//
// use ihash(key) % NReduce to choose the reduce
// task number for each KeyValue emitted by Map.
//
func ihash(key string) int {
	h := fnv.New32a()
	h.Write([]byte(key))
	return int(h.Sum32() & 0x7fffffff)
}

//
// main/mrworker.go calls this function.
//
func Worker(mapf func(string, string) []KeyValue,
	reducef func(string, []string) string) {
	// Your worker implementation here.
	for {
		task := askForTask()
		if task.TaskType == MAP_TYPE {
			mapData(task, mapf)
		} else if task.TaskType == REDUCE_TYPE {
			reduceData(task, reducef)
		} else if task.TaskType == WAIT_TYPE {
			time.Sleep(time.Second)
			continue
		} else {
			return
		}

		call("Master.TaskDone", &task, &TaskCallReply{})
	}
	// uncomment to send the Example RPC to the master.
	//CallExample()
}

func askForTask() TaskCallReply {
	args := TaskCallArgs{}
	reply := TaskCallReply{}
	call("Master.GetTask", &args, &reply)
	return reply
}

func mapData(task TaskCallReply, mapf func(string, string) []KeyValue) {
	filename := task.TaskName
	content := getFileContent(filename)
	kva := mapf(filename, content)
	generateIntermediateFiles(task, kva)
}

func getFileContent(filename string) string {
	file, err := os.Open(filename)
	if err != nil {
		log.Fatalf("cannot open %v", filename)
	}
	content, err := ioutil.ReadAll(file)
	if err != nil {
		log.Fatalf("cannot read %v", filename)
	}
	file.Close()
	return string(content)
}

func generateIntermediateFiles(task TaskCallReply, kva []KeyValue) {
	// kva is divided into buckets , buckets are arrays of key - values
	// for example one entry of one bucket would be : "word" - [1, 1, 1, 1]
	// buckets is Map< int, <Map<string, string[]>> >
	// OR {bucketNum : {key : [values]}}
	buckets := make(map[int]map[string][]string)

	NReduce := task.N
	for i := 0; i < NReduce; i++ {
		buckets[i] = make(map[string][]string)
	}

	for _, kv := range kva {
		key := ihash(kv.Key) % NReduce
		buckets[key][kv.Key] = append(buckets[key][kv.Key], kv.Value)
	}

	for num, kva := range buckets {
		filename := "mr-" + strconv.Itoa(task.TaskNumber) + "-" + strconv.Itoa(num) + ".json"
		dir, _ := os.Getwd()
		file, _ := ioutil.TempFile(dir, filename)
		enc := json.NewEncoder(file)
		enc.Encode(&kva)
		file.Close()
		os.Rename(file.Name(), filename)
	}
}

func reduceData(task TaskCallReply, reducef func(string, []string) string) {
	data := make(map[string][]string) //{key : [values]}

	for X := 0; X < task.N; X++ {
		filename := "mr-" + strconv.Itoa(X) + "-" + task.TaskName + ".json"
		file, _ := os.OpenFile(filename, os.O_CREATE|os.O_RDWR, 0644)
		dec := json.NewDecoder(file)
		var kva map[string][]string //{key : [values]}
		dec.Decode(&kva)
		file.Close()
		for key, values := range kva {
			data[key] = append(data[key], values...)
		}
	}

	outputFile, _ := os.Create("mr-out-" + task.TaskName + ".txt")
	for key, values := range data {
		output := reducef(key, values)
		fmt.Fprintf(outputFile, "%v %v\n", key, output)
	}
	outputFile.Close()
}

//
// example function to show how to make an RPC call to the master.
//
// the RPC argument and reply types are defined in rpc.go.
//
func CallExample() {

	// declare an argument structure.
	args := ExampleArgs{}

	// fill in the argument(s).
	args.X = 99

	// declare a reply structure.
	reply := ExampleReply{}

	// send the RPC request, wait for the reply.
	call("Master.Example", &args, &reply)

	// reply.Y should be 100.
	fmt.Printf("reply.Y %v\n", reply.Y)
}

//
// send an RPC request to the master, wait for the response.
// usually returns true.
// returns false if something goes wrong.
//
func call(rpcname string, args interface{}, reply interface{}) bool {
	// c, err := rpc.DialHTTP("tcp", "127.0.0.1"+":1234")
	sockname := masterSock()
	c, err := rpc.DialHTTP("unix", sockname)
	if err != nil {
		log.Fatal("dialing:", err)
	}
	defer c.Close()

	err = c.Call(rpcname, args, reply)
	if err == nil {
		return true
	}

	fmt.Println(err)
	return false
}
