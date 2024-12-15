import jsonschema
import json
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Validate STRchive loci JSON")
    parser.add_argument("loci", help="STRchive loci JSON file to validate")
    parser.add_argument("schema", help="schema JSON file to validate against")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    # Load the schema
    with open(args.schema) as f:
        schema = json.load(f)

    # Load the data
    with open(args.loci) as f:
        instance = json.load(f)

    # Validate
    any_errors = False
    for locus in instance:
        for field in locus:        
            try:
                jsonschema.validate(locus[field], schema['properties'][field])
            except jsonschema.exceptions.ValidationError as e:
                any_errors = True
                sys.stderr.write(f"JSON data is invalid: {e}\n")
                sys.stderr.write(f"Locus: {locus["id"]}, Field: {field}, Value: {locus[field]}\n")
    if any_errors:
        sys.exit(1)

if __name__ == "__main__":
    main()