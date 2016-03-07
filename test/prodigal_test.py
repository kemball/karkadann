from karkadann.prodigal import *
from karkadann.database import data_location
import unittest as ut 


class TestOverlap(ut.TestCase):
	

	def test_overlap(self):
		from Bio.SeqFeature import FeatureLocation	

		def overlap_exercise(one,two,result):
			self.assertEqual(overlap(one,two),overlap(two,one))
			self.assertEqual(overlap(one,two),result)

		reference_feature= FeatureLocation(100,200,strand=+1)
		inside_feature = FeatureLocation(150,175,strand=+1)
		overlap_exercise(reference_feature,inside_feature,True)

		left_feature = FeatureLocation(75,150,strand=+1)
		overlap_exercise(left_feature,reference_feature,True)

		left_touchy = FeatureLocation(50,100,strand=+1)
		overlap_exercise(left_touchy,reference_feature,False)

		left_inside = FeatureLocation(100,125,strand=+1)
		overlap_exercise(left_inside,reference_feature,True)

		definitely_outside = FeatureLocation(500,1000,strand=-1)
		overlap_exercise(definitely_outside,reference_feature,False)



class FloatingOverlapTest(ut.TestCase):

	def test_left(self):
		from Bio.SeqFeature import FeatureLocation

		def overlap_exercise(one,two,result):
			self.assertEqual(floverlap(one,two),floverlap(two,one))
			self.assertEqual(floverlap(one,two),result)
			if one.strand == 1:
				newstrand = -1
			elif one.strand == 0 or one.strand is None:
				newstrand = 0
			else:
				newstrand = 1
			neg = FeatureLocation(one.start,one.end,newstrand)
			self.assertEqual(floverlap(neg,two),floverlap(one,two))
			self.assertEqual(floverlap(neg,two),result)
		ref = FeatureLocation(100,200,+1)

		left = FeatureLocation(50,75,-1)

		overlap_exercise(ref,left,0)

		on_left = FeatureLocation(50,150,+1)

		overlap_exercise(ref,on_left,.5)

		left_touchy = FeatureLocation(50,100,+1)

		overlap_exercise(left_touchy,ref,0)

		inset_left = FeatureLocation(100,150,+1)

		overlap_exercise(inset_left,ref,.5)

		inside_left = FeatureLocation(125,175,+1)

		overlap_exercise(inside_left,ref,.5)

	def test_right(self):
		from Bio.SeqFeature import FeatureLocation

		def overlap_exercise(one,two,result):
			self.assertEqual(floverlap(one,two),floverlap(two,one))
			self.assertEqual(floverlap(one,two),result)
			if one.strand == 1:
				newstrand = -1
			elif one.strand == 0 or one.strand is None:
				newstrand = 0
			else:
				newstrand = 1
			neg = FeatureLocation(one.start,one.end,newstrand)
			self.assertEqual(floverlap(neg,two),floverlap(one,two))
			self.assertEqual(floverlap(neg,two),result)

		ref = FeatureLocation(100,200,-1)

		right = FeatureLocation(250,300,+1)

		overlap_exercise(ref,right,0)

		on_right = FeatureLocation(200,300,-1)

		overlap_exercise(ref,on_right,0)

	def test_weird(self):
		from Bio.SeqFeature import FeatureLocation

		def overlap_exercise(one,two,result):
			self.assertEqual(floverlap(one,two),floverlap(two,one))
			self.assertEqual(floverlap(one,two),result)
			if one.strand == 1:
				newstrand = -1
			elif one.strand == 0 or one.strand is None:
				newstrand = 0
			else:
				newstrand = 1
			neg = FeatureLocation(one.start,one.end,newstrand)
			self.assertEqual(floverlap(neg,two),floverlap(one,two))
			self.assertEqual(floverlap(neg,two),result)


		ref = FeatureLocation(0,1000)

		odd = FeatureLocation(100,200)

		overlap_exercise(ref,odd,.1)

		very_long = FeatureLocation(500,2000)

		overlap_exercise(very_long,odd,0)

		overlap_exercise(ref,very_long,1/3.0)

		snug = FeatureLocation(1,1000)

		overlap_exercise(ref,snug,999/1000.0)



class TestAnnotate(ut.TestCase):
	#TODO use resources properly, like for the config file. :/
	testgb  = os.path.join(data_location,"test/testassem.gb")
	testrec = SeqIO.parse(testgb,'genbank')
	testrec = list(testrec)
	preserve_testrec = annotate(testrec,preserve_anno=True)
	spoiled_testrec  = annotate(testrec,preserve_anno=False)
	def test_annotate(self):
		#annotation should not delete contigs
		self.assertEqual( len(self.spoiled_testrec) ,len(self.preserve_testrec))
		preserved_contig = self.preserve_testrec[0]
		raw_contig = self.spoiled_testrec[0]
		#or change their order or sequence
		self.assertEqual(len(preserved_contig.seq),len(raw_contig.seq))

	def test_merge(self):
		testcontig = self.testrec[0]
		preservecontig = self.preserve_testrec[0]
		spoiledcontig = self.spoiled_testrec[0]
		testprots = filter(lambda iscds: iscds.type == "CDS",testcontig.features)
		preserveprots = filter(lambda iscds: iscds.type =="CDS",preservecontig.features)
		spoiledprots = filter(lambda iscds: iscds.type =="CDS",spoiledcontig.features)
		self.assertGreaterEqual(len(preserveprots),len(testprots))
		self.assertGreaterEqual(len(preserveprots),len(spoiledprots))
		if testcontig.id == "NZ_KK070022.1":
			print "testprots"
			print testprots
			print "spoiledprots"
			print spoiledprots
			print "preserveprots"
			print preserveprots
			self.assertEqual(len(spoiledprots),3)
			self.assertEqual(len(preserveprots),3)
			self.assertEqual(len(testprots),2)



	

if __name__=="__main__":


	testgb = os.path.join(data_location,"test/testassem.gb")
	testrec=SeqIO.parse(testgb,'genbank')
	betteranno = annotate(testrec)
	SeqIO.write(betteranno,'test.gb','genbank')
	testrec=SeqIO.parse(testgb,'genbank')
	preserved_anno = annotate(testrec,preserve_anno=True)
	SeqIO.write(preserved_anno,'test2.gb','genbank')
	#ut.main()

