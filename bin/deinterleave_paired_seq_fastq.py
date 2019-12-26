import sys

class deinterleave:
    def __init__(self):
        print('Creating an instance of deinterleave class')
        self.interleaved_fastq_path = sys.argv[1]
        print(f'interleaved fastq path is {self.interleaved_fastq_path}')
        self.r1_stream = open(sys.argv[2], 'w')
        self.r2_stream = open(sys.argv[3], 'w')
        print(f'out areguments are {sys.argv[2]} {sys.argv[3]}')

    def deinterleave_fastq(self):
        [self.r1_stream.write(line) if (i % 8 < 4) else self.r2_stream.write(line) for i, line in enumerate(open(self.interleaved_fastq_path, 'r'))]
        self.r1_stream.close()
        self.r2_stream.close()

di = deinterleave()
di.deinterleave_fastq()