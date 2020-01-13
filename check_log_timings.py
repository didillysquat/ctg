import datetime
import calendar
base_names = ['M_18_1922_HUME-ST-35-45_AD002_Lane1','M_18_1922_HUME-ST-35-45_AD002_Lane2','M_18_1923_HUME-ST-5-7_AD007_Lane1', 'M_18_1923_HUME-ST-5-7_AD007_Lane2', 'M_18_1924_HUME-ST-8-11_AD019_Lane1', 'M_18_1924_HUME-ST-8-11_AD019_Lane2', 'M_17_426_1_AD08_Lane1', 'M_17_426_1_AD08_Lane2']
processes_to_check = ['trim_reads', 'rcorrector']
with open('.nextflow.log', 'r') as f:
    lines = [line.rstrip() for line in f]

month_abbrev_to_int =  dict((v,k) for k,v in enumerate(calendar.month_abbr))

def get_time_from_log_line(line):
    month = line.split()[0].split('-')[0]
    if month == 'Dec':
        year = 2019
    else:
        year = 2020
    day = line.split()[0].split('-')[1]
    hour = line.split()[1].split(':')[0]
    minute = line.split()[1].split(':')[1]
    second = line.split()[1].split(':')[2].split('.')[0]
    return datetime.datetime(year=int(year), month=month_abbrev_to_int[month], day=int(day), hour=int(hour), minute=int(minute), second=int(second))

# look through the log file for a line that contains
# basically we want to go through each of the base names and look to see when the processess
# started and when they finished. We want to get a rough idae about how long they took to compelte
# we can have a list of the base names and a list of the processes that we want ot check on
for process in processes_to_check:
    for base_name in base_names:
        # A bool to let us know whether we found an end time or not.
        end_found = False
        start_found = False
        # first look for the start of the process
        for line in lines:
            if 'Submitted process' in line and process in line and base_name in line:
                # then this is the start of the process
                start_time = get_time_from_log_line(line)
                print(f'Start time for {process} for {base_name} is {start_time}')
                start_found = True
            if 'COMPLETED' in line and process in line and base_name in line:
                # then this is when the process was completed
                end_time = get_time_from_log_line(line)

                # Once we have found an end time then we should work out how long this process took
                # and report it
                time_taken = end_time - start_time
                print(f'{process} for {base_name} took {time_taken}\n\n')
                # we don't need to go through the log file anyfurther if we have found this
                end_found = True
                break
        # Here we have either found an end time and reported how long the process took
        # or we have not found an end time.
        # or we have not found a start either.
        # Deal with these possibilites here.
        if not start_found and not end_found:
            print(f'{process} for {base_name} has not started yet\n\n')
        elif start_found and not end_found:
            print(f'{process} for {base_name} has been running for {datetime.datetime.now() - start_time} so far and is still not complete\n\n')