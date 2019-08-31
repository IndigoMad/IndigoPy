#!/usr/bin/python3
import IndigoPy.Protein.pdb as pdb

a=pdb.Protein('4hqp')
insert=a[0]['J'][25:43]
b=pdb.Protein('1WN9')
c=b.couple('A',35,46,insert,name='1WN9Insert')
c.exportpdb()
