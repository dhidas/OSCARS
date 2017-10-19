#!/bin/bash

jupyter nbconvert --to python Example_*ipynb GettingStarted_*ipynb --TemplateExporter.exclude_markdown=True --TemplateExporter.exclude_input_prompt=True

FILES=`ls Example_*py | grep -v OSG`
for i in $FILES; do sed -i '1,6d' $i; done

