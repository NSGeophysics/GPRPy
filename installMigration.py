from subprocess import call
import os
import shutil

print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('You are downloading an external software with a more restrictive license\n')
print('If you are a student or researcher using it for non-commercial purposes,\nthen you should be fine. Otherwise, please see the license here:\n')
print('https://www.crewes.org/ResearchLinks/FreeSoftware/\n')
print('and here:\n')
print('https://github.com/njwilson23/irlib\n')

try:
    call(['git','clone','https://github.com/AlainPlattner/irlib.git'])
    shutil.move('irlib','gprpy/irlib')
except:
    print('no git installed, simply downloading')
    os.mkdir('irlib')
    os.mkdir('irlib/external')
    import urllib.request
    urllib.request.urlretrieve ("https://raw.githubusercontent.com/AlainPlattner/irlib/master/external/mig_fk.py", "irlib/external/mig_fk.py")
    urllib.request.urlretrieve ("https://raw.githubusercontent.com/AlainPlattner/irlib/master/__init__.py", "irlib/__init__.py")
    shutil.move('irlib','gprpy/irlib')
