from timeit import default_timer as timer
import parameters_start

def main():
    start = timer()
    parameters_start.main()
    end = timer()
    print('Time:',end-start)

if __name__ == "__main__":
    main()