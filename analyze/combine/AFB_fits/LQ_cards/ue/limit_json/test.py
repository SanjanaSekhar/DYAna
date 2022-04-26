import json

with open('limits_eu.json', 'r+') as f:
    data = json.load(f)
    data['1500.0']['exp+1'] = 123 # <--- add `id` value.
    #f.seek(0)        # <--- should reset file position to the beginning.
    json.dump(data, f,indent=3)
    f.truncate()     # remove remaining part
