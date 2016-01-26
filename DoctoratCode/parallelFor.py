from joblib import Parallel, delayed
import multiprocessing

def processInputs(i):
    return i*i


if __name__=='__main__':
    inputs=range(10)
    print inputs
    num_cores=multiprocessing.cpu_count()
    print num_cores
    result=Parallel(num_cores)(delayed(processInputs)(i) for i in inputs)

