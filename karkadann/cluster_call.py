from database import Hit,Gene,Cluster
from collections import defaultdict

def _classify(gene):
	s = defaultdict(int)
	for score, hmm in gene.hit_scores():
		s[hmm] = score
	PKS_KS = s["PKS_KS"]
	if PKS_KS>50 and PKS_KS > s["bt1fas"] and PKS_KS > s["ft1fas"] and PKS_KS > s["hg1D"] and PKS_KS > s['hg1E'] and PKS_KS > s['fabH']:
		return "PKS_I"
	if s['tra_KS']>65 and s['PKS_AT'] < 50 < PKS_KS:
		return "PKS_I_transAT"
	if s['t2ks']>50:
		if max([s['t2ks'], s['ene_KS'], s['mod_KS'], s['hyb_KS'], s['itr_KS'], s['tra_KS'], s['bt1fas'], s['ft1fas'], s['t2fas'], s['hg1d'], s['hg1e'], s['fabH']]) == s['t2ks']:
			return "PKS_II"
	if s['t2clf']>50:
		if max([s['t2clf'], s['ene_KS'], s['mod_KS'], s['hyb_KS'], s['itr_KS'], s['tra_KS'], s['bt1fas'], s['ft1fas'], s['t2fas'], s['hg1D'], s['hg1E'], s['fabH']])==s['t2clf']:
			return "PKS_II"
	if s['Chal_sti_synt_C']>35 or s['Chal_sti_synt_N']>35:
		return "PKS_III"
	if s['Condensation']>20 and (s['AMP_binding']>20 or s['A_OX']>20):
		return 'nrps'
	if s['Terpene_synth_C']>23 or s['phytoene_synt']>20 or s['Lycopene_cycl']>80 or s['terpene_cyclase']>50 or s['NapT7']>250 or s['fung_ggpps']>420 or s['fung_gpps2']>312 or s['dmat'] or s['trichodiene_synth']>150:
		return 'terpene'
	if s['LANC_like']>80 or s['Land_dehyd_N']>20 or s['Lant_dehyd_C']>20 or s['TIGR03731']>18:
		return 'lanti'
	if s['BLS']>250 or s['CAS']>250 or s['tabtoxin']>500:
		return 'beta_lactam'
	if s['strH_like']>50 or s['strK_like1']>800 or s['strK_like2']>650 or s['neoL_like']>50 or s['DOIS']>500 or s['valA_like']>600 or s['spcFG_like']>200 or s['spcDK_like_cou']>600 or s['salQ']>480:
		return 'amino'
	if s['IucA_IucC']>30:
		return 'sidero'
	if s['ectoine_synt']>35:
		return 'ectoine'
	if s['AfsA']>25:
		return 'butryo'
	if s['indsynth']>100:
		return 'indole'
	if s['LipM']>50 or s['LipU']>30 or s['LipV']>375 or s['ToyB']>175 or s['TunD']>200 or s['pur6']>200 or s['pur10']>600 or s['nikJ']>200 or s['nikO']>400:
		return 'nucleo'
	if s['MoeO5']>65:
		return "phosphoglycolipids"
	if s['melC']>40:
		return "melanin"
	if s['ppm']>150:
		return 'phos'
	if s['nos']>150:
		return "nos"
	if s['lan_c_jd']>25:
		return 'lant_c_jd'
	if s['terpene_jd']>40:
		return 'terpene_jd'
	if s['TOMM']>25:
		return "TOMMdocking"


def call_clusters(contig):
	genes = list(contig.genes())
	classifications = list(map(_classify, genes))
	clusters = defaultdict(list)
	for i, gc in enumerate(zip(genes, classifications)):
		g, c = gc
		if c:
			clusters[c].append(genes[i-6:i+7])
	for kind in clusters.keys():
		for i, listed in enumerate(sorted(clusters[kind], key=lambda x: x[0].location.start)):
			if not len(listed) or listed == clusters[kind][-1]:
				continue
			if len(clusters[kind])==1: #no self-overlap
				continue
			if listed[-1] in clusters[kind][i+1]:
				# we have overlap. do something?
				clusters[kind][i]=list(set(clusters[kind][i]) | set(clusters[kind[i+1]]))
				clusters[kind][i+1] = []
	final_clusters = []
	for kind in clusters.keys():
		if len(clusters[kind]):
			for genes in clusters[kind]:
				new_cluster = Cluster(gene_list=genes, classification=kind)
				new_cluster.save()
				final_clusters.append(new_cluster)
	return final_clusters




