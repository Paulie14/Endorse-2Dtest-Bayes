from flow_wrapper import Wrapper

if __name__ == "__main__":
    wrap = Wrapper(solver_id=1)
    wrap.set_parameters(data_par=[6e-15, 0.17])
    res = wrap.get_observations()

