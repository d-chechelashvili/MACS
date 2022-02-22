from subprocess import Popen, PIPE, run

server1 = Popen(['python3', './project_3_2/rsa_signing_server.py', '49104'])
server2 = Popen(['python3', './project_3_2/rsa_verify_server.py', '49105'])

p = run(['python3', 'project_3_2/sign.py'], stdout=PIPE, stderr=PIPE)
print(p.stdout)

server1.terminate()
server2.terminate()
