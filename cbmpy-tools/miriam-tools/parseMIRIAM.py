import os, time, numpy, re, pprint
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
#import pyscescbm as cbm

PP = pprint.PrettyPrinter()

## mDir = os.path.join(cDir, 'miriam')
mDir = cDir
mFile = os.path.join(mDir,'miriamresources.xml')

F = file(mFile,'r')
res = F.read()
F.close()

datatype = re.compile('<datatype id=.*?</datatype>', re.DOTALL)

res2 = re.findall(datatype, res)
print len(res2)


name = re.compile('<name>.*?</name>')
uri = re.compile('<uri type="URL">.*?</uri>')
dataEntry = re.compile('<dataEntry>.*?</dataEntry>', re.DOTALL)
regex = re.compile('pattern=.*?>')
exam = re.compile('<dataEntityExample>.*?</dataEntityExample>')

out = {}
for r_ in res2:
    nm = re.findall(name, r_)[0].replace('<name>','').replace('</name>','').strip()
    ur = re.findall(uri, r_)[0].replace('<uri type="URL">','').replace('</uri>','').strip()
    de = re.findall(dataEntry, r_)[0].replace('<dataEntry>','').replace('</dataEntry>','').strip()
    rx = re.findall(regex, r_)[0].replace('pattern="','').replace('">','').strip()
    ex = re.findall(exam, r_)
    if len(ex) > 0:
        ex = ex[0].replace('<dataEntityExample>','').replace('</dataEntityExample>','').strip()
    else:
        ex = None
    out.update({ nm : {
        'name' : nm,
        'url' : ur,
        'data_entry' : de,
        'pattern' : rx,
        'example' : ex
        }})

p = PP.pformat(out)
print p

F = file(os.path.join(mDir, 'miriamids.py'), 'w')
F.write('# created on %s\n\n' % time.strftime('%y%m%d:%H%M'))
F.write('miriamids =\\\n')
F.write(p)
F.write('\n')
F.close()
