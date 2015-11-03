import os

version = raw_input()

s1 = os.system("git tag %s")
if s1 == 0:
   s2 = os.system("git tag -f latest")

if s1 == 0 and s2 == 0:
   s2 = os.system("git push --tags -f")
   

