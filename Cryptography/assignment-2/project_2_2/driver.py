from subprocess import Popen, PIPE, run

server2 = Popen(['./project_2_2/mac-server', '6667'])
server3 = Popen(['./project_2_2/vrfy-server', '6668'])

p = run(['python3', 'project_2_2/forge.py', 'project_2_2/challenge_message.txt'], stdout=PIPE, stderr=PIPE, encoding='utf-8')
print(p.stdout)

server2.terminate()
server3.terminate()
