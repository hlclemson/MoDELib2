import re

# Example text
text = """
periodicDipoleNodes=10 10; # number of extra nodes on each dipole
periodicDipolePoints=-150.0 0.0 0.0
                      -50.0 0.0 0.0;
periodicDipoleHeights=2550 2550; # height of each dipole, in number of planes
"""

# New coordinates to replace the old ones
new_coordinates = "per = -100.0 0.0 0.0\n-200.0 0.0 0.0"

# Regular expression to find and replace the line with 'periodicDipolePoints'
#updated_text = re.sub(r'(periodicDipolePoints=)[-\d.\s]+;', r'\1' + new_coordinates + ';', text)
updated_text = re.sub(r'periodicDipolePoints=[-\d.\s]+;',  new_coordinates, text)

print(updated_text)
