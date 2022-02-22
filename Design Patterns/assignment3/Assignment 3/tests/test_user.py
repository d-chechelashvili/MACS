from app.user import User


def test_create() -> None:
    username = "user1"
    user = User(username)
    assert user is not None


def test_get_name() -> None:
    username = "user1"
    user = User(username)
    assert user.get_name() == username


#
# def test_on_video_published() -> None:
#     username = "user1"
#     user = User(username)
#     captured_output = io.StringIO()
#     sys.stdout = captured_output
#     user.on_video_published()
#     sys.stdout = sys.__stdout__
#     assert captured_output.getvalue() == username + "\n"
