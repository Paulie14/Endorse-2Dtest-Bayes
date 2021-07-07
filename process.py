from flow_wrapper import Wrapper
import time

if __name__ == "__main__":
    wrap = Wrapper(solver_id=1)
    wrap.set_parameters(data_par=[6e-15, 0.17])
    t = time.time()
    res = wrap.get_observations()
    print(res)
    print("LEN:", len(res[1]))
    print("TIME:", time.time()-t)
