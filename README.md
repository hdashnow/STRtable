# STRchive

("ess-tee-archive")

Short Tandem Repeat disease loci resource

[⭐️ View the data at strchive.org ⭐️](http://strchive.org/)

<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="http://strchive.org/">STRchive</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://github.com/hdashnow">Harriet Dashnow</a> is licensed under <a href="http://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1"></a></p>

## Contributors

- Harriet Dashnow
- Laurel Hiatt
- Akshay Avvaru
- Vincent Rubinetti

## Contributing

If you notice an error, omission, or update, feel free to leave a comment or create a pull request.

To make a change to the STRchive data itself, please edit `data/STRchive-loci.json`

Then run the "linting" script and fix any errors:  
`python scripts/check-loci.py data/STRchive-loci.json`

## Development

### Update TRGT genotyping catalogs

```
python scripts/make-catalog.py -g hg38 -f TRGT data/STRchive-loci.json data/hg38.STRchive-disease-loci.TRGT.bed
python scripts/make-catalog.py -g T2T -f TRGT data/STRchive-loci.json data/T2T-chm13.STRchive-disease-loci.TRGT.bed
python scripts/make-catalog.py -g hg19 -f TRGT data/STRchive-loci.json data/hg19.STRchive-disease-loci.TRGT.bed
```

### Update extended BED files

```
python scripts/make-catalog.py -f bed -g hg38 data/STRchive-loci.json data/hg38.STRchive-disease-loci.bed
python scripts/make-catalog.py -f bed -g T2T data/STRchive-loci.json data/T2T-chm13.STRchive-disease-loci.bed
python scripts/make-catalog.py -f bed -g hg19 data/STRchive-loci.json data/hg19.STRchive-disease-loci.bed
```