from timeit import default_timer as timer
import parameters_start

start = timer()
parameters_start.main()
end = timer()
print('Time:',end-start)