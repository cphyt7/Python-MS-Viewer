import xml.etree.ElementTree as ET
import gzip
import sys
import csv
from base64 import b64decode
from array import array
import matplotlib.pyplot as plt

#Inputs and file handling
try:
    spectrafilename = sys.argv[1]
except NameError:
    print >>sys.stderr, 'Spectra file not found. Please include spectra file as argument 1.'
    sys.exit(1)

try:
    inputmasstable = open(sys.argv[2])
except NameError:
    print >>sys.stderr, 'Mass table not found. Please include mass table as argument 2.'
    sys.exit(1)
except IOError:
    print >>sys.stderr, 'Mass table file not found. Please ensure file name spelled correctly.'
    sys.exit(1)

try:
    scan = sys.argv[3]
except IndexError:
    print >>sys.stderr, 'Scan number not found. Please include scan number as argument 3.'
    sys.exit(1)

try:
    inputsequence = sys.argv[4].upper()
except IndexError:
    print >>sys.stderr, 'Protein sequence not found. Please include protein sequence as argument 4.'
    sys.exit(1)

masstable = csv.reader(inputmasstable)

try:
    seqfile = gzip.open(spectrafilename)
except IOError:
    print >>sys.stderr, 'Spectra file not found. Please ensure file name spelled correctly.'
    sys.exit(1)
    
#Creating dictionary to store mass information
massdict = {}
for row in masstable:
    key = row[0]
    massdict[key] = row[1]

#Parsing for correct scan
document = ET.iterparse(seqfile)
ns = "{http://sashimi.sourceforge.net/schema/}"

for event,ele in document:
    if ele.tag == ns+"scan":
        if ele.get('num') == scan:
            peakselt = ele.find(ns+'peaks')
            peaks = array('f',b64decode(peakselt.text))
            if sys.byteorder != 'big':
                peaks.byteswap()
            mzs = peaks[::2]
            ints = peaks[1::2]

#Computing b-ions and y-ions from input sequence
reversedsequence = inputsequence[::-1]
ionset = set()
print
print'***B-ions***\n'

btotalmass = 1
bmass = 0
bcount = 0
totalions = {}
totalsequences = {}
try:
    for amino in inputsequence:
        bmass = float(massdict[amino])
        btotalmass += bmass
        bcount +=1
        key = btotalmass
        totalions[key] = 'b'+str(bcount)
        key2 = 'b'+str(bcount)
        totalsequences[key2] = inputsequence[:bcount]
        ionset.add(btotalmass)
        print 'The mass of',inputsequence[:bcount],'b'+str(bcount),'is:',btotalmass
except KeyError:
    print >>sys.stderr, 'Peptide sequence contains incorrect amino acids. Please re-enter sequence.'
    sys.exit(1)
print
print '***Y-ions***\n'

ytotalmass = 19
ymass = 0
ycount = 0
for amino in reversedsequence:
    ymass = float(massdict[amino])
    ytotalmass += ymass
    ycount +=1
    key = ytotalmass
    totalions[key] = 'y'+str(ycount)
    key2 = 'y'+str(ycount)
    totalsequences[key2]= reversedsequence[:ycount]
    ionset.add(ytotalmass)
    print 'The mass of',reversedsequence[:ycount],'y'+str(ycount),'is:',ytotalmass

#Matching, figure generation and annotations
mzint = {}

try:
    normalizedints = map(lambda x: (x/max(ints)*100),ints)
except NameError:
    print >>sys.stderr, 'Scan number not found. Please try another scan number.'
    sys.exit(1)
    
mzint = dict(zip(mzs,normalizedints))

figure = plt.stem(mzs,normalizedints,markerfmt=' ')
plt.title("MS/MS for scan number: "+scan+', Matching peptide: '+inputsequence)
count2 = 0
count3 = 0
for k, v in mzint.items():
    for j in ionset:
        if k-0.3 < j < k+0.3 and v > 5.0:
            if totalions[j][0] == 'y':
                plt.annotate(totalions[j]+'\n'+str(j),xy=(k,v),ha ='center',color='maroon')
                plt.figtext(.905,.85-count2,(totalions[j]+': '+totalsequences[totalions[j]]),fontdict=None)
                count2 += 0.05
            else:
                plt.annotate(totalions[j]+'\n'+str(j),xy=(k,v),ha ='center',color='navy')
                plt.figtext(.905,.85-count2,(totalions.get(j)+': '+totalsequences[totalions[j]]),fontdict=None)
                count2 += 0.05


plt.xlabel('M/Z')
plt.ylabel('Intensity %')
plt.ylim(0,110)
plt.figtext(.905,.9,'Ion fragment sequences:')
plt.setp(figure,color='black')

plt.show()


