from Bio import Entrez

with open('/home/maycon/Documents/LAB/TS/data/raw/vert_sequence_id', 'r') as ids_file:
    ids = [each_id.strip() for each_id in ids_file.readlines()]

    nuc_ids = []
    missing_nuc = []

    for each_id in ids:
        Entrez.email = 'flayner5@gmail.com'
        search = Entrez.esearch(db='nuccore', term=each_id)
        result = Entrez.read(search)
        print(f'{each_id} | {result}')
        if result['IdList']:
            nuc_id = result['IdList'][0]
            nuc_ids.append((each_id, nuc_id))
        else:
            missing_nuc.append(each_id)

    out_file = open('/home/maycon/Documents/LAB/TS/data/vert_taxids.txt', 'w')
    out_file.write('accession,taxid\n')

    for each_acc, each_id in nuc_ids:
        Entrez.email = 'flayner5@gmail.com'
        fetch = Entrez.efetch(db='nuccore', id=each_id, retmode='xml')
        fetch_res = Entrez.read(fetch)
        print(fetch_res)

        for dic in fetch_res[0]['GBSeq_feature-table'][0]['GBFeature_quals']:
            if 'taxon:' in dic['GBQualifier_value']:
                out_file.write(
                    f"{each_acc},{dic['GBQualifier_value'].replace('taxon:', '')}\n")

    out_file.close()

    print(f'Missing nuc id: {missing_nuc}')
