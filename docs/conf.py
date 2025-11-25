import sys
from pathlib import Path

# Add the '_ext' directory to sys.path
# sys.path.append(str(Path('_ext').resolve()))
# sys.path.insert(0, str(Path(__file__).parent / "_ext"))

sys.path.append(str(Path('_ext').resolve()))


# Import the echoautobody extension
# extensions = ['echoautobody', 'helloworld']
extensions = ['echoautobody']

# Manually import extension (this registers the directives)
# import echoautobody
# import helloworld