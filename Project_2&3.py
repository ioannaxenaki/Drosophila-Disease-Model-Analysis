### Drosophila Disease Model Analysis Project (Project 2 & 3)

############# Start Project 2 #############
 
import pandas as pd

flybase_file="/home/ioanna/metaptixiako/kanterakis/disease_model_annotations_fb_2025_02.tsv"
mondo_file="/home/ioanna/metaptixiako/kanterakis/mondo.json"



####Task 1
keep_cols=[0,3,4,5]
col_names=['FBgn ID','DO qualifier','DO ID','DO term']

data=pd.read_csv(flybase_file, sep='\t', comment='#', usecols=keep_cols, names=col_names, dtype=str) #data frame with 4 columns (col_names)

data.dropna(inplace=True)

data_rows=len(data)
unique_genes=len(data['FBgn ID'].unique())
unique_qualifiers=data['DO qualifier'].unique()
do_ids=set(data['DO ID'])
unique_diseases=len(data['DO term'].unique())



####Task 2
import json

#dictionary to match do ids to mondo
do_ids_to_mondo={}

with open(mondo_file, 'r') as f:
    mondo_data=json.load(f)

nodes=mondo_data['graphs'][0]['nodes'] #list with nodes

for node in nodes:
    node_id=node.get('id')
    if not node_id or not node_id.startswith("http://purl.obolibrary.org/obo/MONDO_"): #we want only nodes with MONDO
        continue
    
    meta=node.get('meta') #we check if the node has meta
    if meta:
        xrefs=meta.get('xrefs') #we check if meta has xrefs
        if xrefs:
            for xref in xrefs:
                xref_val=xref.get('val') #we get the val of each xref
                if xref_val and xref_val.startswith('DOID:'): #we only need vals with DOID
                    if xref_val not in do_ids_to_mondo: #if we see the DOID (key) for the first time we make an empty list (value) 
                        do_ids_to_mondo[xref_val]=[]

                    do_ids_to_mondo[xref_val].append(node_id) #each DOID (key) has a list of all the vals (value)

print(f'We found {len(do_ids_to_mondo)} matches DOID in mondo.json')

#we want to get only MONDO IDs that we have the matching DO ID from flybase
mondo_ids=set()
doids_matched_count=0
doids_unmatched_count=0

for doid in do_ids: #unique do ids from flybase
    if doid in do_ids_to_mondo: #we check if our DOID is in the dict with {DOID:MONDO}
        mondo_ids.update(do_ids_to_mondo[doid])
        doids_matched_count+=1
    else:
        doids_unmatched_count+=1

print(f'Unique DO IDs from flybase file: {len(do_ids)}')
print(f'DO IDs that matched to MONDO: {doids_matched_count}')
print(f'DO IDs that did not match to MONDO: {doids_unmatched_count}')
print(f'Unique MONDO Ids: {len(mondo_ids)}')



####Task 3
#to do task 3 we need the functions of the project 3

############# Pause Project 2 Task 3 #############

############# Start Project 3 #############

###function part_1
def part_1(json_file_path): #reads mondo.json and creates dictionary that holds the tree
    tree={'labels':{}, 'category':{}, 'disease':{}} 

    with open(json_file_path, 'r') as f:
        data=json.load(f)
    
    nodes=data['graphs'][0]['nodes'] #each node is a disease
    edges=data['graphs'][0]['edges'] #they show the relation between the connected nodes. each edge desribe a relation 

    for node in nodes: #we get the labels for each node
        node_id=node.get('id')
        if node_id and node_id.startswith("http://purl.obolibrary.org/obo/MONDO_"):
            if 'lbl' in node:
                tree['labels'][node_id]=node['lbl']

    
    for edge in edges:
        if edge.get('pred')=="is_a": #'pred' is the kind of relation. 'is_a' means 'is a kind of' 
            disease_id=edge.get('sub') #'sub' makes the relation (Leukemia)
            category_id=edge.get('obj') #'obj' accepts the relation (Cancer)

            if (disease_id and disease_id.startswith("http://purl.obolibrary.org/obo/MONDO_") and category_id and category_id.startswith("http://purl.obolibrary.org/obo/MONDO_")):
                #for each disease
                if disease_id not in tree['category']:
                    tree['category'][disease_id]=[]
                tree['category'][disease_id].append(category_id) #for each disease we have all categories
                #for each category
                if category_id not in tree['disease']:
                    tree['disease'][category_id]=[]
                tree['disease'][category_id].append(disease_id) #for each category we have all diseases

    return tree

###function part_2
def part_2(mondo_tree, mondo_id): #finds labels in the 3rd label of the tree of a MONDO ID
    level_1="http://purl.obolibrary.org/obo/MONDO_0000001"
    #we want to find all nodes in level 3
    if 'level_3_nodes' not in mondo_tree:
        level_3_nodes=set()
        #level 2
        level_2_nodes=mondo_tree['disease'].get(level_1,[])
        #level 3
        for level_2_node in level_2_nodes:
            diseases_of_level_2=mondo_tree['disease'].get(level_2_node,[])
            level_3_nodes.update(diseases_of_level_2) #we have all level 3 nodes
        
        mondo_tree['level_3_nodes']=level_3_nodes
        #print(f'We have {len(level_3_nodes)} nodes in level 3') #we have 98 nodes in level 3

    #we want to find all the above levels/disriptions ('parents') of the given ID disease
    given_id=f"http://purl.obolibrary.org/obo/{mondo_id}"
    all_parents=set() #all the parents of the given id from all the levels

    parents_to_check=[given_id] #we want to check all the ids that are above the given. we save them in the list to check them all until none remains
    parents_checked={given_id} #we save all the discriptions/'parents' that we already found are related to the given id

    while len(parents_to_check)>0: #we check the ids in the list until its empty
        current_id=parents_to_check.pop(0)  #we begin to check with the first id
        parents=mondo_tree['category'].get(current_id,[]) #we get the parents of the first id
        for parent_id in parents: #we check if we have seen each parent
            if parent_id not in parents_checked: #if not
                parents_checked.add(parent_id) #we save it as checked(already seen)
                all_parents.add(parent_id) #we save it as one of all parents of the given id
                parents_to_check.append(parent_id) #we save it to check it later for its own parents
    
    level_3_ids=all_parents.intersection(mondo_tree['level_3_nodes']) #we find the common ids between all the parents of the given id and the ids that are in level 3 of the tree

    #we have to get the label of these label 3 ids
    level_3_labels=[]
    for level_3_id in level_3_ids: #for each id in level 3
        if level_3_id in mondo_tree['labels']: #we check if the id has a label in the mondo tree
            level_3_labels.append(mondo_tree['labels'][level_3_id]) #we save the label of the 3rd level
        
    return list(set(level_3_labels)) #we get the unique labels of the 3rd level about the given id

mondo_tree=part_1(mondo_file)
mondo_categories=part_2(mondo_tree, 'MONDO_0011719')
print(f'The labels in the 3rd level of MONDO_0011719 are: {mondo_categories}') #given example

############# End Project 3 #############

############# Continue Project 2 Task 3 #############

gene_to_categories={} #dicctionary that matches gene from flybase to categories from mondo tree
already_searched_mondo={} #we save searched mondo ids so we dont call part 2 again
all_categories=set() #we save all the categories that we found

unique_gene_with_do_ids=data[['FBgn ID','DO ID']].drop_duplicates() #from the data (datafram from flybase) we get the FBgn ID and the DOID

for number,row in unique_gene_with_do_ids.iterrows(): #for each gene
    gene_id=row['FBgn ID']
    do_id=row['DO ID']

    mondo_id_list=do_ids_to_mondo.get(do_id,[]) #we get the MONDO IDs for this DOID
    gene_categories=set() #categories for this gene

    for mondo_full_id in mondo_id_list:
        mondo_id=mondo_full_id.split('/')[-1] #from the full https address of the MONDO we keep only the MONDO_X

        if mondo_id not in already_searched_mondo: #we check if we have seen this mondo id before
            categories=part_2(mondo_tree, mondo_id) #if not, we find its categories using part_2
            already_searched_mondo[mondo_id]=categories #we save the categories of each mondo id 

        categories=already_searched_mondo[mondo_id] #we wget the categories of each mondo id

        gene_categories.update(categories) #we save the categories found for this gene
        all_categories.update(categories) #we save all the categories we found

    #we update the final dictionary that matches gene to categories
    if len(gene_categories)>0:
        if gene_id not in gene_to_categories: 
            gene_to_categories[gene_id]=set()
        gene_to_categories[gene_id].update(gene_categories)

print(f'We found {len(all_categories)} unique categories from level 3')
print(f'We matched {len(gene_to_categories)} genes to categories')



####Task 4

print(f'We have {unique_genes} genes in total')

#we want to find which genes match to each DO qualifier 
genes_per_qualifier={}

for qualifier in unique_qualifiers: #for each qualifier
    genes=set(data[data['DO qualifier']==qualifier]['FBgn ID']) #we found the genes that match this qualifier
    genes_per_qualifier[qualifier]=genes #we save the genes for each qualifier

#we want to find which genes match to each mondo category
genes_per_category={}

for gene,categories in gene_to_categories.items(): #for each gene and its categories 
    for category in categories: #for each of its categories
        if category not in genes_per_category:
            genes_per_category[category]=set()
        genes_per_category[category].add(gene) #we find which genes belong in this category
        

contingency_tables=[] #list of contigency tables
for category in all_categories: #for each mondo category
    for qualifier in unique_qualifiers: #for each qualifier
        genes_with_category=genes_per_category.get(category,set()) #we get genes in this category
        genes_with_qualifier=genes_per_qualifier.get(qualifier,set()) #we get genes with this qualifier

        b=len(genes_with_category.intersection(genes_with_qualifier)) #genes that have both the category and the qualifier
        c=len(genes_with_category-genes_with_qualifier) #genes that have ony the category
        q=len(genes_with_qualifier-genes_with_category) #genes that have only the qualifier
        n=unique_genes-(b+c+q) #genes that have neither the category or the qualifier

        contingency_table=[[b,c],[q,n]] #the contigency table
        contingency_tables.append({'category':category, 'qualifier':qualifier, 'table':contingency_table, 'b':b, 'c':c, 'q':q, 'n':n})

print(f'We have {len(contingency_tables)} contigency tables')


####Task 5
from scipy.stats import fisher_exact

for res in contingency_tables: #for each result (its a dictionary) in the list of tables
    table=res['table'] #we get the table
    odds_ratio, p_value=fisher_exact(table, alternative='two-sided') #we get the odds ratio and p-value for each table by computing the two-sided test
    res['odds_ratio']=odds_ratio #how possible to see the relation between the category and the qualifier (OR=1:a gene with this category is equally possible to have the qualifier as a gene without the category (there is no relation), OR>1:a gene with this category is more likely to have this qualifier than a gene without the category (positive relation), OR<1:a gene with this cetegory is less likely to heve this qualifier than a gene without the category (negative relation))
    res['p_value']=p_value #the possibility to get numbers equally or more extreme that the results we got if the null hypothesis (there is no relation) is true


####Task 6
contingency_tables.sort(key=lambda x: x['p_value']) #we sort by the p-value
top_10_results=contingency_tables[:10] #we get the 10 lowest results

for res in top_10_results:
    p_val=res['p_value']
    odd_ra=res['odds_ratio']
    categ=res['category']
    qual=res['qualifier']
    print(f'p-value:{p_val},    odds ratio:{odd_ra},   category:{categ},    qualifier:{qual}')


####Bonus 1
#the expected value of genes that have both the category and the qualifier is: ((b+c)*(b+q))/unique_genes
#the fold_change is: b/expected value

for res in contingency_tables:
    b=res['b']
    c=res['c']
    q=res['q']
    total_category=b+c
    total_qualifier=b+q

    expected_b=(total_category*total_qualifier)/unique_genes

    if expected_b>0:
        fold_change=b/expected_b
    else:
        fold_change=0

    res['expected_b']=expected_b
    res['fold_change']=fold_change

contingency_tables.sort(key=lambda x: x['p_value']) #we sort by the p-value
top_10_results=contingency_tables[:10] #we get the 10 lowest results

for res in top_10_results:
    p_val=res['p_value']
    odd_ra=res['odds_ratio']
    categ=res['category']
    qual=res['qualifier']
    exp_val=res['expected_b']
    fold_ch=res['fold_change']
    observed=res['b']
    print(f'p-value:{p_val},    odds ratio:{odd_ra},   category:{categ},    qualifier:{qual},     observed b:{observed},     expected b:{exp_val},      fold change:{fold_ch}')


###Bonus 2
from statsmodels.stats.multitest import multipletests

p_values=[res['p_value'] for res in contingency_tables] #list of p-values

reject, pvals_corrected, _, _=multipletests(p_values, alpha=0.05, method='fdr_bh')
#reject is a list of True/False, if it still rejects the null hypothesis
#pvals_corrected is a list with the adjusted p-values

for i in range(len(contingency_tables)): #we add for each table
    contingency_tables[i]['p_adjusted']=pvals_corrected[i] #the adjusted p-value
    contingency_tables[i]['rejects_null']=reject[i] #if it still rejects the null hypothesis

count=0
for res in contingency_tables:
    if res['rejects_null']:
        p_val=res['p_value']
        odd_ra=res['odds_ratio']
        categ=res['category']
        qual=res['qualifier']
        p_adj=res['p_adjusted']
        fold_ch=res['fold_change']
        print(f'p-value:{p_val},    odds rario:{odd_ra},   category:{categ},    qualifier:{qual},     adjusted p-value:{p_adj},    fold change:{fold_ch}')
        count+=1
print(f'There are {count} results that still reject the null hypothesis after multiple test correction')