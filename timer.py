from timeit import default_timer as timer
import parameters_start
import probable_throwshed

def main():
    start = timer()
    parameters_start.main()
    # probable_throwshed.main()
    end = timer()
    print('Time:',end-start)

if __name__ == "__main__":
    main()