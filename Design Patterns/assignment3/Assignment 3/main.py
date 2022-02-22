from app.repl import CLI, ICLI
from app.repository import FileDBRepository
from app.video_sharing_system import SubscriberPrinter, VideoSharingSystem


def simulate_repl(repl: ICLI) -> None:
    print("Enter empty command to break")
    while True:
        user_input = input(">>> ")
        if user_input == "":
            break
        command = repl.parse_cmd(user_input)
        repl.run(command)


if __name__ == "__main__":
    system = VideoSharingSystem(FileDBRepository())
    system.attach(SubscriberPrinter())
    repl = CLI(system)
    simulate_repl(repl)
