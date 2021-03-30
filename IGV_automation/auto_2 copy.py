import socket

from pysam import VariantFile

vcf_in = VariantFile(b"/Users/smgreen02/Documents/broad/klebs/results/vcf/frags_final_variants.vcf")  # auto-detect input format
vcf_out = VariantFile('-', 'w', header=vcf_in.header)

HOST = 'localhost'
PORT = 60151

#socket = IGVsocket(host='localhost', port=60151)
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((HOST, PORT))

s.sendall(b'new\n')

from_server = s.recv(4096)

print(from_server)

s.sendall(b'genome /Users/smgreen02/Documents/broad/klebs/klebs_new/reference/Klebsiella_pneumoniae_subsp_pneumoniae_MGH_78578_plasmids.fasta\n')

from_server = s.recv(4096)

print(from_server)

s.sendall(b'load /Users/smgreen02/Documents/broad/klebs/results/bam/nanopore.aligned.sorted.bam\n')

from_server = s.recv(4096)

print(from_server)

s.sendall(b'load /Users/smgreen02/Documents/broad/klebs/results/vcf/frags_final_variants.vcf\n')

from_server = s.recv(4096)

print(from_server)

for rec in vcf_in.fetch():
    print(rec.contig)
    s.sendall(b'goto ' + str.encode(rec.contig) + b'\n')
    from_server = s.recv(4096)
    print(from_server)
    s.sendall(b'snapshotDirectory /Users/smgreen02/Documents/broad/klebs/results/bam/snapshotDirec/screenshots\n')
    from_server = s.recv(4096)
    print(from_server)
    s.sendall(b'snapshot ' + str.encode(rec.contig) + b'\n')
    from_server = s.recv(4096)
    print(from_server)

s.sendall(b'exit\n')

from_server = s.recv(4096)

print(from_server)

s.close()