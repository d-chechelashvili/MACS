import sys
import getopt
from collections import defaultdict

import Checksum
import BasicSender

'''
This is a skeleton sender class. Create a fantastic transport protocol here.
'''


class Sender(BasicSender.BasicSender):
    def __init__(self, dest, port, filename, debug=False, sackMode=False):
        super(Sender, self).__init__(dest, port, filename, debug)
        self.sackMode = sackMode
        self.debug = debug
        self.packet_size = 1400
        self.window_size = 7
        self.fast_retransmit_limit = 4
        self.timeout = 0.5
        self.ack_counter = defaultdict(int)
        self.truly_received = -1  # denotes index
        self.seqno = 0

    # Main sending loop.
    def start(self):
        self.stop_and_wait_send("syn")
        self.transmit_data()
        self.stop_and_wait_send("fin")
        self.infile.close()

    def stop_and_wait_send(self, msg_type):
        packet = self.make_packet(msg_type, self.seqno, "")
        while True:
            self.send(packet)
            ack = self.receive(timeout=self.timeout)
            if self.correct_ack(ack): break
        self.seqno += 1
        self.truly_received += 1

    def correct_ack(self, ack):
        if not ack: return False
        _, seqno, _, _ = self.split_packet(ack)
        seqno = int(seqno.split(';', 1)[0])
        return seqno == (self.seqno + 1) and Checksum.validate_checksum(ack)

    def transmit_data(self):
        num_to_read = self.window_size
        num_to_send = num_to_read
        packets = []
        received_packets = []
        while True:
            new_packets, has_more_data = self.read_data(num_to_read)
            packets += new_packets
            if self.sackMode: received_packets += [False for _ in new_packets]
            if not has_more_data and self.truly_received == len(packets): break
            self.send_packets(packets, num_to_send, received_packets)
            ack = self.receive(timeout=self.timeout)
            num_to_read, num_to_send = self.slide_window(ack, received_packets)

    def read_data(self, num):
        packets = []
        has_more_data = True
        for _ in xrange(num):
            data = self.infile.read(self.packet_size)
            if not data:
                has_more_data = False
                break
            packet = self.make_packet("dat", self.seqno, data)
            self.seqno += 1
            packets.append(packet)
        return packets, has_more_data

    def send_packets(self, packets, num_to_send, received_packets):
        packets_to_send = packets[self.truly_received: self.truly_received + num_to_send]
        for i, packet in enumerate(packets_to_send):
            if not self.sackMode or not received_packets[self.truly_received + i]:
                self.send(packet)

    def slide_window(self, ack, received_packets):
        num_to_read = 0
        if ack is None:
            num_to_send = self.window_size
        else:
            _, raw_seqno, _, _ = self.split_packet(ack)
            seqno = int(raw_seqno.split(';', 1)[0])
            if seqno <= self.truly_received:
                num_to_send = 0
            elif not Checksum.validate_checksum(ack):
                num_to_send = 1 #0
            else:
                if self.sackMode:
                    received_packets_indices = [
                        int(packet) for packet in raw_seqno.split(';', 1)[1].split(',') if packet
                    ]
                    for index in received_packets_indices:
                        received_packets[index - 1] = True
                self.ack_counter[seqno] += 1
                if self.ack_counter[seqno] == self.fast_retransmit_limit:
                    num_to_send = 1
                    self.ack_counter[seqno] = 0
                    # num_to_send = self.window_size
                    # for i in xrange(self.truly_received + 1, self.truly_received + self.window_size + 1):
                    #     self.ack_counter[i] = 0
                else:
                    num_to_read = seqno - 1 - self.truly_received
                    num_to_send = num_to_read #self.window_size
                    self.truly_received = seqno - 1
        return num_to_read, num_to_send

'''
This will be run if you run this script from the command line. You should not
change any of this; the grader may rely on the behavior here to test your
submission.
'''
if __name__ == "__main__":
    def usage():
        print "BEARS-TP Sender"
        print "-f FILE | --file=FILE The file to transfer; if empty reads from STDIN"
        print "-p PORT | --port=PORT The destination port, defaults to 33122"
        print "-a ADDRESS | --address=ADDRESS The receiver address or hostname, defaults to localhost"
        print "-d | --debug Print debug messages"
        print "-h | --help Print this usage message"
        print "-k | --sack Enable selective acknowledgement mode"


    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   "f:p:a:dk", ["file=", "port=", "address=", "debug=", "sack="])
    except:
        usage()
        exit()

    port = 33122
    dest = "localhost"
    filename = None
    debug = False
    sackMode = False

    for o, a in opts:
        if o in ("-f", "--file="):
            filename = a
        elif o in ("-p", "--port="):
            port = int(a)
        elif o in ("-a", "--address="):
            dest = a
        elif o in ("-d", "--debug="):
            debug = True
        elif o in ("-k", "--sack="):
            sackMode = True

    s = Sender(dest, port, filename, debug, sackMode)
    try:
        s.start()
    except (KeyboardInterrupt, SystemExit):
        exit()
