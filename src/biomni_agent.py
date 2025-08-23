import sys
sys.path.append("../")

from biomni.agent import A1


agent= A1(path="data/bundle/SlAREB1_bundle_SGN_first/", llm="gpt-oss:120b")

logs = agent.go("I have CDS of tomato(solanum lycopersicum) AREB gene locus located in data/bundle/SlAREB1_bundle_SGN_first/CDS", \
"I have cDNA of tomato(solanum lycopersicum) AREB gene locus located in data/bundle/SlAREB1_bundle_SGN_first/cDNA", \
"I have the genomic sequence of tomato(solanum lycopersicum) AREB gene locus located in data/bundle/SlAREB1_bundle_SGN_first/genomic", \
"I hope to design high-quality CRISPR guide RNA for gene editing and subsequent wet-lab validation in order to confer saline tolerance in tomato", \
"I have the enumerated guides and reinforcement learning-selected stored in data/bundle/SlAREB1_bundle_SGN_first/guides", \
"I have stored some preliminary results from the experimental setup in data/bundle/SlAREB1_bundle_SGN_first/figs/", \
"I have additional genomic data on tomato (solanum lycopersicum) located in data/", \
"I have a summary of project implementation located in data/bundle/project_summary",\
"Please help me with the next steps in communicating results and validating findings in my project with the aim to confer saline tolerance in gene-edited tomato using CRISPR gene editing. I have all the data and resources stored in the specified paths.")
