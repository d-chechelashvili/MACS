from subprocess import Popen, PIPE, run

server1 = Popen(['./project_2_1/server', '6667'])

p = run(['python3', 'project_2_1/decipher.py', 'project_2_1/ctext.txt'], stdout=PIPE, stderr=PIPE)
print(p.stdout)

server1.terminate()
