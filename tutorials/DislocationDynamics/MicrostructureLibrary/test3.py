import re

# Example text
#text = """
#periodicDipoleNodes=10 10; # number of extra nodes on each dipole
#periodicDipolePoints=120 0 0
#                    -100.0 0.0 0.0;
#periodicDipoleHeights=2550 2550; # height of each dipole, in number of planes
#"""
text = """
periodicDipolePoints=120 0 0
                    -100.0 0.0 0.0;
"""

# Corrected regular expression substitution
updated_text = re.sub(r'periodicDipolePoints=((?:.|\s)*?);', 'periodicDipolePoints=100;', text)

print(updated_text)
