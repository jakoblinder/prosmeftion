import os
import fcntl
import subprocess
import time
import re

from pathlib import Path
from typing import Optional


class PyFORM:

    pipes: list
    process: subprocess.Popen
    prompt: str
    n_bytes: int

    def __init__(
        self,
        form_path: Path,
        n_pipes: int,
        prompt: str = "",
        n_bytes: int = 512,
        input_dir: Path = None,
        *args,
    ):
        """
        Establishes a pipe connection to FORM for a fast, two-way communication of Python and FORM without the generation of input and output files.
        According to the 'Embedding FORM in other applications' chapter from the FORM manual (https://www.nikhef.nl/~form/maindir/documentation/reference/online/online.html) two pipes need to be generated,
            r_form, w_py = os.pipe()  # Py -> FORM
            r_py, w_form = os.pipe()  # FORM -> Py
        with one 'going from' Python to FORM and one from FORM to Python.
        The file descriptors (fds) 'r_form' and 'w_form' are used to read messages send from Python into the fd 'w_py' and send messages to Python which can there be read with the fd 'r_py'.
        In order establish the connection, FORM needs to know its fds and, therefore, has to be started with the command line option
            -pipe r_form,w_form
        . In order to check the correctness of the connection, some process ids are send back and forth.
        Of course, more than one communication channel can be made, as specified by the number of pipe pairs, 'n_pipes'.

        ----------
        form_path
            Path of the FORM program to be run.
        n_pipes
            Number of pipes the FORM requires.
        prompt
            Prompt used in FORM to determine the end of an expression (make sure thats true!!!).
        n_bytes
            Maximum number of bytes that can be read in one step from a FORM message.
        input_dir
            FORM is executed with the "-p" option and a path of a directory for input, include, procedure and subroutine files.
        args (optional)
            Additional options FORM should be run with.

        Returns
        -------
        """
        self.prompt = prompt
        self.n_bytes = n_bytes
        self.pipes = [os.pipe() + os.pipe() for i in range(n_pipes)]  # Create pairs of pipes.
        pipe_ids = ",".join(f"{r_py},{w_py}" for r_py, _, _, w_py in self.pipes)
        # "Removes the .str file on crash, whatever its contents. (Under ordinary circumstances at a crash a .str file will not be removed if it has a nonzero content.)"
        # command = ["form", "-pipe", pipe_ids]
        command = ["form"]
        command.append("-Z")
        if input_dir:
            command.append("-p")
            command.append(str(input_dir))
        command.extend(args)
        command.append(form_path.name)

        # Try to set environment variable:
        os.environ["FORM_PIPES"] = pipe_ids

        self.process = subprocess.Popen(
            command,
            pass_fds=[fd for pipe in self.pipes for fd in pipe],
            cwd=form_path.resolve().parent,
            # stdin=subprocess.DEVNULL,  # .PIPE
            stdout=subprocess.PIPE,  # .PIPE
            stderr=subprocess.STDOUT,  # .STDOUT
            # FIXME: Handle somehow external channel error. Maybe catch the error and just start FORM again?
        )
        # Check correctness of connection by sending process ids back and forth.
        for _, w_form, r_form, _ in self.pipes:
            assert (
                int(os.read(r_form, self.n_bytes).decode().strip()) == self.process.pid
            )
            os.write(
                w_form, f"{self.process.pid},{os.getpid()}\n".encode(encoding="ascii")
            )

        # for stdout_line in iter(self.process.stdout.readline, ""):
        #     print(stdout_line.decode())
        # time.sleep(1e-3)
        # with self.process.stdout as output:
        #     ou = output.read(5).decode()
        #     print(ou)
        #     if re.search(r"Error initializing preset external channels", ou):
        #         print("AHAAA")

        time.sleep(1e-10)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        # Wait for process to be finished (otherwise FORM terminates within execution).
        self.return_code = self.process.wait()
        # while self.process.poll() is None:
        #     time.sleep(1e-5)
        self.close()
        if self.return_code:
            raise ConnectionError("FORM didn't initiated the pipes successfully.")

    def close(self):
        """
        End connection and close all pipes.
        """
        self.process.terminate()
        for fds in self.pipes:
            for fd in fds:
                try:
                    os.close(fd)
                except OSError:
                    pass

    def write(self, n_pipe: int, content: str, prompt: Optional[str] = None):
        """
        Send message to FORM.
        The message is end with one linebreak '\n' followed by the prompt as specified by the '#prompt' preprocessor variable and then an additional linebreak '\n'.
            '\n<prompt>\n'
        The default prompt is an empty string.
        The prompt can be changed with the command
            #prompt <prompt>
        which cam be specified differently for every pipe.
        -------
        n_pipe
            Specify in which pipe the message is send.
        content
            The message itself.
        prompt
            A different prompt.
        Returns
        -------
        """
        fd = self.pipes[n_pipe - 1][1]
        os.write(fd, f"{content}\n{prompt or self.prompt}\n".encode(encoding="ascii"))

    def read(self, n_pipe: int, n_bytes: Optional[int] = None):
        """
        Read message send from FORM into pipe 'n_pipe'.
        -------
        n_pipe
            Specify from which pipe the message should be read.
        n_bytes
            Maximum number of bytes of the message.
        -------
        Returns
            Message send from FORM in the given pipe.
        """
        # Wait for process to be finished
        self.process.wait()
        # while self.process.poll() is None:
        #     time.sleep(1e-5)

        fd = self.pipes[n_pipe - 1][2]
        return os.read(fd, n_bytes or self.n_bytes).decode()

    def readline(self, n_pipe: int):
        """
        Read one line from the message send from FORM into pipe 'n_pipe'.
        Requires that the message ends with an '\n', since FORM doesn't append a linebreak by itself.
        Otherwise, it will deadlock.
        -------
        n_pipe
            Specify from which pipe the message should be read.
        -------
        Returns
            On line of the message send from FORM in the given pipe.
        """
        # Wait for process to be finished
        self.process.wait()
        # while self.process.poll() is None:
        #     time.sleep(1e-5)

        fd = self.pipes[n_pipe - 1][2]
        with os.fdopen(fd, closefd=False) as readFile:
            return readFile.readline()

    def read_all(self, n_pipe: int):
        """
        Read everything send from FORM into pipe 'n_pipe'.
        -------
        n_pipe
            Specify from which pipe the message should be read.
        -------
        Returns
            All messages send from FORM in the given pipe.
        """
        # Wait for process to be finished
        self.process.wait()
        # while self.process.poll() is None:
        #     time.sleep(1e-5)

        fd = self.pipes[n_pipe - 1][2]
        # only UNIX: https://stackoverflow.com/questions/31113805/reading-from-unbuffered-os-fdopen-file-object-doesnt-behave-like-os-read
        fcntl.fcntl(fd, fcntl.F_SETFL, os.O_NONBLOCK)
        with os.fdopen(fd, 'rb', buffering=0, closefd=False) as readFile:
            content = readFile.read()
            if content:
                return content.decode()
            else:
                return content

if __name__ == "__main__":
    form_path = Path("pyform.frm")
    with PyFORM(form_path, 2) as form:
        form.write(1, "Symbol a;\n\nSet something: a;\n", prompt="READY")
        form.write(1, "3*a*a", prompt="READY")
        form.write(2, "3*b")

        res1 = form.readline(1)
        res2 = form.readline(2)

    print("1: ", res1)
    print("2: ", res2)
