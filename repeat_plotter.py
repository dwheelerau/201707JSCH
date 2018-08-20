from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import seaborn as sns
import sys

from ConfigParser import SafeConfigParser


# SETTINGS
# NOTE 1.
# make this a sys.argv[1]
tar_fasta = sys.argv[1] #'./CST20R1_24_0.fasta' # './testfile_2.fasta' 

# config file
configFilePath = "config.ini"
config = SafeConfigParser()
config.read(configFilePath)

#print config.get('main', 'key1') # -> "value1"
#print config.get('main', 'key2') # -> "value2"
#print config.get('main', 'key3') # -> "value3"

AA_unit = int(config.get('main', 'AA_unit')) #2 
dna_unit = AA_unit * 3
AA_skipNterm = int(config.get('main','AA_skipNterm')) #8 
dna_skipNterm = AA_skipNterm * 3
# If set AA_skipCterm to 0 if you want to capture all of the seq at the 3' end
AA_skipCterm = int(config.get('main', 'AA_skipCterm')) #0 # 12:-24 = 24 OR set to 0 to include all 3' seq
dna_skipCterm = AA_skipCterm * 3 # will be set to 0 if above is 0

strains = {}

with open(tar_fasta) as f:
    for rec in SeqIO.parse(f, 'fasta'):
        # DNA skip from AA skip
        data = []
        dna_start = rec[:dna_skipNterm]
        # logic to capture entire seq if requried
        if dna_skipCterm == 0:
            dna_string = rec[dna_skipNterm:] # captures 3' end
        else:
            dna_string = rec[dna_skipNterm:-dna_skipCterm]   
        data.append(('START', 'START'))
        start = 0
        end = dna_unit
        while end <= len(dna_string):
            dna_chunk = dna_string[start:end]
            # table 12 is yeast alt
            data.append((str(dna_chunk.seq), str(dna_chunk.translate(table=12).seq)))
            start = end
            end = end + dna_unit
        data.append(('END', 'END' ))
        strains[rec.id] = data
        
print(strains)
# NOTE 3. This needs to be saved as an outfile somehow

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 10)) #9,4

# colect all repeats and find unique ones to assign a color
uq_repeats = list(set([repeat[1] for repeat in strains for repeat in strains[repeat]]))

#np.random.seed(sum(map(ord, "palettes")))

col_dict = {}
col_opts = sns.xkcd_rgb.keys()
n = 0
for rep in uq_repeats:
    if rep == 'XXXXX':
        # set bad repeat to black, redundent
        col_dict[rep] = (0,0,0)
    else:
        col_dict[rep] = sns.xkcd_rgb[np.random.choice(col_opts, replace=False)] #palet[n]
    n+=1

## repeat plot
start_end = []
width = 1
spacer = 0

for i,strain in enumerate(strains.keys()):
    used_colors_dict = {}
    used_colors = []
    start = 0
    for dna, repeat in strains[strain]:
        if start == 0:
            # start of seq (not repeat)
            col_dict[repeat] = '#000000'
            col = '#000000'
            start_end.append(repeat)        

        elif repeat == strains[strain][-1][1]:
            # end of seq (not repeat)
            col_dict[repeat] = '#000000'
            col = '#000000'
            start_end.append(repeat)
        else:
            col=col_dict[repeat]
            spacer = 0
        start += width
        if col not in used_colors or col == '#ffffff':
            used_colors.append(col)
            used_colors_dict[col] = [(dna, repeat)]
            hatch=""
        else:
            # pep same, but DNA diff, make hatch
            found,notfound = 0,0
            for d, r in used_colors_dict[col]:
                try:
                    assert r == repeat
                except AssertionError:
                    print('start, end or error')
                    print(r)
                    print(repeat)
                if dna != d:
                    notfound+=1
                else:
                    found+=1
            used_colors_dict[col].append((dna, repeat))
            # ocnce it is found, dont hatch it again
            #if notfound != 0 and found ==0:
            # jan seems to hash even if it has been found before
            if notfound != 0:
                hatch="/"
            else:
                hatch = ""
        # could put logic here to make it white box
        axes[0].barh(i,width, left=start - width, color=col_dict[repeat],
                     hatch=hatch, align='center', edgecolor='black')

# get rid of xtics
axes[0].set_xticklabels('')
ind = np.arange(len(strains))
axes[0].set_yticks(ind)
labels1 = strains.keys()
axes[0].set_yticklabels(labels1, fontsize=8) # SETS Font size

# remove ticks and borders
axes[0].get_xaxis().set_ticks([])
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[0].spines['bottom'].set_visible(False)
axes[0].spines['left'].set_visible(False)
axes[0].tick_params(axis=u'both', which=u'both',length=0)

# remove the non-repeat start and ends before making the key
assert len(start_end) % 2 == 0
for norep in start_end:
    col_dict.pop(norep, None)
    
# custom key on right plot
col2 = col_dict.values()
labels2 = col_dict.keys()
axes[1].barh(np.arange(len(labels2)),len(labels2)*[0.3], left=0.1, color=col2, align='center')
axes[1].set_xbound((0,1))
axes[1].set_xticklabels('')

# fix ticks
axes[1].set_yticks(np.arange(len(labels2)))
axes[1].set_yticklabels(labels2)
axes[1].get_xaxis().set_ticks([])
#plt.tick_params(top='off', bottom='off', left='on', right='off', labelleft='on', labelbottom='on')
plt.tick_params(top=False, bottom=False, left=True, right=False, labelleft=True, labelbottom=True)

# turn of border and ticks
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].spines['bottom'].set_visible(False)
axes[1].spines['left'].set_visible(False)
axes[1].tick_params(axis=u'both', which=u'both',length=0)

# # NOTE 2. set outfile as sys.argv
fig.tight_layout()
fig.savefig('outfile.pdf')
plt.show()
