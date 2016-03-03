from karkadann.prodigal import *
from karkadann.database import data_location

def overlap_test():
	from Bio.SeqFeature import FeatureLocation

	def test_overlap(one,two,result):
		assert(overlap(one,two)==overlap(two,one))
		assert(overlap(one,two)==result)

	reference_feature= FeatureLocation(100,200,strand=+1)
	inside_feature = FeatureLocation(150,175,strand=+1)
	test_overlap(reference_feature,inside_feature,True)

	left_feature = FeatureLocation(75,150,strand=+1)
	test_overlap(left_feature,reference_feature,True)

	left_touchy = FeatureLocation(50,100,strand=+1)
	test_overlap(left_touchy,reference_feature,False)

	left_inside = FeatureLocation(100,125,strand=+1)
	test_overlap(left_inside,reference_feature,True)

	definitely_outside = FeatureLocation(500,1000,strand=-1)
	test_overlap(definitely_outside,reference_feature,False)

if __name__=="__main__":
	#I should really convert this to use the .cfg file
	# why else did I write it?
	overlap_test()
	testgb = os.path.join(data_location,"test/testassem.gb")
	testrec=SeqIO.parse(testgb,'genbank')
	betteranno = annotate(testrec)
	SeqIO.write(betteranno,'test.gb','genbank')
	testrec=SeqIO.parse(testgb,'genbank')
	preserved_anno = annotate(testrec,preserve_anno=True)
	SeqIO.write(preserved_anno,'test2.gb','genbank')
