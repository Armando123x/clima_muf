gfortran -o iri iritest.for irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for rocdrift.for

python /scripts/refresh_indices.py
python /scripts/update_images.py
 
sleep infinity