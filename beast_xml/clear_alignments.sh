#!/bin/bash

for file in *.xml; do
    python3 -c "
import re, sys

with open('$file', 'r') as f:
    content = f.read()

cleaned = re.sub(r'<alignment\b.*?</alignment>', '', content, flags=re.DOTALL)

with open('$file', 'w') as f:
    f.write(cleaned)
" 
    echo "Processed: $file"
done