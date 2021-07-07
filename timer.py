import time


class TimeTimer:
    def __init__(self, check_with_print=False):
        self.created_time = time.time()
        self.checked_time_list = []
        self.time_call_counter = 0
        self.check_with_print = check_with_print
        self.checked_time_list.append(self.created_time)

    def check_time(self, msg=""):
        current_time = time.time()
        new_time = current_time
        last_time = self.checked_time_list[-1]

        self.time_call_counter += 1
        self.checked_time_list.append(new_time)

        if self.check_with_print:
            interval = new_time - last_time
            info = f"time info -> count: [{self.time_call_counter}] elapsed time: [{interval}]"
            if len(msg) > 0:
                info = f"time info -> count: [{self.time_call_counter}] msg: {msg} elapsed time: [{interval}]"
            print(info)

