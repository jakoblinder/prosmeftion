import os
import re
import subprocess
import time
from pathlib import Path

from typing import Optional, List


def pyForm(form_path: Path, write: Optional[str] = [], prompt: Optional[str] = "", input_dir: Optional[Path] = None, preprocessor_only=False, debug=False, recursion_depth=0) -> str:
    """
    Establishes a pipe connection to FORM for a fast, two-way communication of Python and FORM without the
    generation of input and output files.
    According to the 'Embedding FORM in other applications' chapter from the FORM manual
    (https://www.nikhef.nl/~form/maindir/documentation/reference/online/online.html) two pipes need to be generated,
        r_form, w_py = os.pipe()  # Py -> FORM
        r_py, w_form = os.pipe()  # FORM -> Py
    with one 'going from' Python to FORM and one from FORM to Python.
    The file descriptors (fds) 'r_form' and 'w_form' are used to read messages send from Python into the fd 'w_py'
    and send messages to Python which can there be read with the fd 'r_py'. In order establish the connection,
    FORM needs to know its fds and, therefore, has to be started with the command line option
        -pipe r_form,w_form
    . In order to check the correctness of the connection, some process ids are sent back and forth.

    In order to avoid deadlocks caused by FORM, all messages FORM should receive are sent in the beginning and FORMs
    response is read only once in the end.

    Write in FORM: "#setexternal `PIPE1_'" to specify the Pipe.
    Messages are received in form with "#fromexternal"
    and written to Python with e.g. "#toexternal "%E\n", expression".

    :param form_path:
        Path of the FORM program to be run.
    :param write:
        List of Expressions which are piped to FORM.
    :param prompt:
        Prompt (e.g.: "READY") used in FORM (-> "#prompt READY")to determine the end of an expression.

    :param input_dir:
        FORM is executed with the "-p" option and a path of a directory for input, include, procedure
        and subroutine files.

    :return:
        Messages send by FORM into the pipe.
    """
    if debug:
        with open(form_path) as file:
            formprogram = file.read()
            formprogram = re.sub(r"#setexternal `PIPE1_'", "", formprogram)
            formprogram = re.sub(r"#prompt READY", "", formprogram)
            formprogram = re.sub(r"#toexternal", "#write", formprogram)
            for message in write:
                formprogram = re.sub(r"#fromexternal",message,formprogram, count=1)
        return formprogram

    r_py, w_form = os.pipe()
    r_form, w_py = os.pipe()
    extra_args = ["-q", "-M"]  # ["-q", "-M"]
    if input_dir:
        extra_args += ["-p", str(input_dir)]
    if preprocessor_only:
        # Run only the preprocessor and dump its output.
        extra_args.append("-y")
    command = (
        ["form", "-pipe", f"{r_py},{w_py}"]
        + extra_args
        + [form_path]
    )
    # Establish connection to FORM:
    with subprocess.Popen(
        command,
        pass_fds=(r_py, w_form, r_form, w_py),
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    ) as process:
        os.close(r_py)
        os.close(w_py)
        os.write(w_form, f"{process.pid},{os.getpid()}\n".encode(encoding="ascii"))
        assert int(os.read(r_form, 8).decode()) == process.pid
        try:
            for i, message in enumerate(write):
                os.write(w_form, f"{message}\n{prompt}\n".encode(encoding="ascii"))
                if preprocessor_only:
                    print(f"{i+1}. Message:")
                    print(f"{message}\n{prompt}\n".encode(encoding="ascii"))
        except BrokenPipeError:
            return pyForm(form_path, write, prompt, input_dir)
        finally:
            os.close(w_form)

        if process.stdout:
            console_output = process.stdout.read().decode()
            if re.search(r"terminating", console_output) or preprocessor_only:
                # FORM has thrown an error
                raise subprocess.SubprocessError(f"FORM message:\n{console_output}")

    with os.fdopen(r_form, "rb") as f:
        data = f.read1()
        if data:
            while data1 := f.read1():
                data += data1
            data = data.decode()
            data = re.sub(r"\\", "", data)  # Remove newline character from FORM.
            return data
        else:
            if recursion_depth >= 10:
                print(f"Aborted after {recursion_depth + 1} tries.")
                raise subprocess.SubprocessError(f"FORM message:\n{console_output}")
            print("New_Try!")
            return pyForm(form_path, write, prompt, input_dir, recursion_depth=recursion_depth+1)

if __name__ == "__main__":
    form_path = Path("pyform_1channel.frm").resolve()
    expression = 4095*"1" + "5"  # "3*a*a"
    write = ["Symbol a;\n\nSet something: a;\n",
             expression]
    for i in range(10000):
        print(i + 1)
        res1 = pyForm(form_path, write, "READY")

    print("1: ", res1)