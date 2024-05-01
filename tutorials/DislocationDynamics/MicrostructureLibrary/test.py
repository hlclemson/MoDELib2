import re

# Path to your text file
file_path = 'periodicDipole.txt'

# Open and read the file
with open(file_path, 'r') as file:
    text = file.read()

# Regular expression to find the line with 'periodicDipolePoints'
match = re.search(r'periodicDipolePoints=([-\d.\s]+);', text)

print("Match:", match.group(1).strip())
exit()
# Check if a match was found
if match:
    # Extract the coordinates and split into lines
    coordinates_block = match.group(1).strip()
    # Split the block into individual coordinate strings
    coordinates_list = coordinates_block.split('\n')
    # Strip extra spaces and split each line into a list of floats
    coordinates = [list(map(float, coords.strip().split())) for coords in coordinates_list if coords.strip()]
    print("Coordinates for periodicDipolePoints:", coordinates)
else:
    print("No 'periodicDipolePoints' found in the file.")
