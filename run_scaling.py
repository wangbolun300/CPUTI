import subprocess
import os
import sys

def scaling(n_threads):
    path = os.getcwd()
    exe = os.path.join(path, "test_bin")

    if sys.platform == "win32":
        exe += ".exe"

    cmd = [exe] + [str(n_threads)]


    out = subprocess.check_output(cmd)
    out = out.decode("utf-8")

    index = out.find("timing temp")
    index+=len("timing temp")
    end = out.find("\n", index)
    algo_time = float(out[index:end])
    end+=2

    index = out.find("timing merge")
    index+=len("timing merge")
    end = out.find("size", index)
    time_merge = float(out[index:end])
    end+=2

    index = out.find("timing sort")
    index+=len("timing sort")
    end = out.find("\n", index)
    sort_time = float(out[index:end])
    
    return algo_time, time_merge, sort_time


if __name__ == "__main__":
    n_trials = 10
    n_threads = [1,2,4,8,16,32]
    out = []
    for t in n_threads:
        tot_algo_time = 0
        tot_time_merge = 0
        tot_sort_time = 0

        print("running {} thread".format(t))
        
        for run in range(n_trials):
            print("\t run {}/{}".format(run+1, n_trials))
            algo_time, time_merge, sort_time = scaling(t)

            tot_algo_time += algo_time
            tot_time_merge += time_merge
            tot_sort_time += sort_time
        
        tot_algo_time /= n_trials
        tot_time_merge /= n_trials
        tot_sort_time /= n_trials

        out.append((t, tot_algo_time, tot_time_merge, tot_sort_time))

    print(out)