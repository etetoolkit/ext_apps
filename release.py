import os

version = raw_input("Version number?")
version = version.strip()
print "version", version
s1 = os.system("git tag %s" %version)
if s1 == 0:
   s2 = os.system("git tag -f latest")

if s1 == 0 and s2 == 0:
   s2 = os.system("git push --tags -f")
   

