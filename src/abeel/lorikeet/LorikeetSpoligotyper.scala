package abeel.lorikeet

import net.sf.samtools.SAMRecordIterator
import net.sf.samtools.BAMFileReader

import scala.collection.JavaConversions._
import net.sf.samtools.SAMFileReader
import java.io.File

import be.abeel.util.CountMap

import java.io.PrintWriter

import org.arabidopsis.ahocorasick.AhoCorasick
import org.arabidopsis.ahocorasick.SearchResult
import net.sf.jannot.utils.SequenceTools
import net.sf.jannot.refseq.Sequence
import net.sf.jannot.refseq.MemorySequence

import scala.io.Source
import atk.util.Tool
import org.apache.commons.math3.distribution.PoissonDistribution
import java.text.NumberFormat
import java.util.Locale

import org.apache.commons.math3.distribution.ExponentialDistribution

import scala.util.Random
import atk.compbio.fastq.FastQFile

import scala.collection.AbstractIterator
import atk.util.TimeInterval

/**
 * Features:
 * - 100% matches only
 * - 43	 spacers
 * - runs entire BAM file
 *
 */
object LorikeetSpoligotyper extends Tool {

  def revcomp(read: Array[Byte]) = {
    val out = Array.ofDim[Byte](read.length)
    for (i <- 0 until read.length) {
      out(read.length - i - 1) = SequenceTools.complement(read(i).toChar).toByte
    }
    out
  }
  case class Config(val spacerFile: String = null, val outputFile: String = null, files: List[File] = List(), val debug: Boolean = false)

  def expand(s: String) = {
    val alpha = List('A', 'T', 'G', 'C')
    val expanded = for {
      i <- 0 until s.length()
      val arr = s.toArray
      val list = for {
        a <- alpha

      } yield {
        arr(i) = a;
        new String(arr)
      }
    } yield (list)
    expanded.flatten.toSet.toList
  }

  /**
   * args(0) = output file
   *
   *
   */
  def main(args: Array[String]): Unit = {

    val genomeLen = 4411708;
    val alignPercent = 0.9
    val pValueThreshold = 0.01
    val parser = new scopt.OptionParser[Config]("java -jar lorikeet.jar spoligotype") {
      opt[String]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("Required: File where you want the output to be written")
      opt[String]('s', "spacer") action { (x, c) => c.copy(spacerFile = x) } text ("Optional: File containing spacers.") 
      opt[Unit]("debug") action { (x, c) => c.copy(debug = true) } text ("Optional: Show debug output.")
      arg[File]("<file>...") unbounded () required () action { (x, c) => c.copy(files = c.files :+ x) } text ("Required: Input files. BAM, SAM, fastq and fastq.gz format are supported.")

    }
    parser.parse(args, Config()) map { config =>
      /* Load spacers */
      val lines = if (config.spacerFile != null) tLines(config.spacerFile).toList else scala.io.Source.fromInputStream(LorikeetSpoligotyper.getClass().getResourceAsStream("/longSpacers.fasta")).getLines().toList; //listlines("v:/TB-ARC/references/spoligoTypeSpacers.fasta")
      val pw = if (config.outputFile != null) new PrintWriter(config.outputFile) else new PrintWriter(System.out)

      pw.println(generatorInfo)

      pw.println("# Alternative spacer input: " + config.spacerFile)
      pw.println("## If the above value is null, the default 43 spacers are used")
      val in = (lines.grouped(2).map(f => (f(0).split("_")(1), f(1))).toList)

      val repeatSequences = List(("repeat", "GTTTCCGTCCCCTCTCGGGGTTTTGGGTCTGACGA"), ("left_repeat", "GTTTCCGTCCCC"), ("middle_repeat", "TCTCGGGGTTTT"), ("right_repeat", "GGGTCTGACGA"))

      
      /**
       * Shuffle spacers as random markers
       */
      val rg = new Random(107)
      val shuffled = in.zipWithIndex.map(triplet => {

        "shuffle_" + triplet._2 -> new String(rg.shuffle(triplet._1._2.toArray.toList).toArray)

      })

      val forwardSpacersLimit = in ++ repeatSequences // ++ mirus()

      /** Expand spacers with all 1 mismatch possibilities */

      val forwardSpacers = forwardSpacersLimit.map(pair => {
        val many = expand(pair._2)
        val zipped = many.zipWithIndex
        zipped.map(m => pair._1 -> m._1)
      }).flatten ++ shuffled

      val rcSpacers = forwardSpacers.map { case (x, y) => (x + "_rc", new String(revcomp(y.getBytes))) }

      val spacers1 = forwardSpacers ++ rcSpacers

      val spacers2 = spacers1.map(f => f._1.replaceAll("d", "") -> f._2)

      val all = spacers2.map(_._2)
     
      val unique = all.groupBy(identity).mapValues(_.size)

      val spacers = spacers2.filter(f => unique(f._2) == 1)
    
      pw.println("# Input files: " + config.files)

      /* Prep index */
      val tree = new AhoCorasick();
      for ((id, spacer) <- spacers) {
        tree.add(spacer.getBytes(), id)
      }
      tree.prepare();

      val nf = NumberFormat.getInstance(Locale.US)
      nf.setMaximumFractionDigits(6)

      for (inputFile <- config.files) {
        pw.println("# Processing: " + inputFile)
        val it = SequenceIterator(inputFile)

        val cm = new CountMap[String]
        /* Adds pseudocounts to all markers */
        val labels = spacers.map(_._1).toSet
        for (sp <- labels) {
          cm.count(sp)
        }

        var progress: Int = 0
        var totalCoverage: Long = 0

        val time = System.currentTimeMillis();
        val nano = System.nanoTime()
        if (config.debug)
          pw.println("# Progress\tExpected coverage\tspeed\ttime\tcurrent progress")
        while (it.hasNext()) {
          val sr = it.next()

          val read = sr.seq

          totalCoverage += sr.seq.length

          val result = tree.search(read.getBytes)

          for (a <- result) {
            for (s <- a.asInstanceOf[SearchResult].getOutputs()) {
              cm.count(s.asInstanceOf[String])
            }
            //       
          }

          // 2 snapshot coverage for all spacers through-out bam file
          progress += 1
          if (progress % 100000 == 0) {
            val expectedCoverage = (totalCoverage * alignPercent) / genomeLen
            print(".")
            if (config.debug)
              pw.println("# " + progress + "\t" + expectedCoverage + "\t" + ((System.nanoTime() - nano) / progress / 1.0e6) + " ms/read\t" + new TimeInterval(System.currentTimeMillis() - time) + "\t" + cm.take(3))

          }

        }

        /* close bam file */
//        it.close
        /* report results */
        if (config.debug)
          pw.println("# Final count map: " + progress + "\t" + cm)
        val expectedCoverage = (totalCoverage * alignPercent) / genomeLen
        /* Take the shuffled tags */
        val tags = cm.filter(p => p._1.matches("shuffle_.*"))
        val list = tags.values.map(_.toInt).toList.sortBy(identity)
        if (config.debug)
          pw.println("List size: " + list.size)
     
        val avg = if (list.size > 0) list.sum.toDouble / list.size else 1
        if (config.debug)
          pw.println("# Sum = " + list.sum)
        if (config.debug)
          pw.println("# average = " + avg)
        val exp = new ExponentialDistribution(avg + 0.1)

        val buffer = new StringBuffer()
        pw.println("# Marker\tread-depth\tp-value\tA/P")
         for (i <- 1 to 43) {
          val z = cm.get("" + i) + cm.get(i + "_rc")
          val p = (1 - exp.cumulativeProbability(z)) * 43
          pw.println(i + "\t" + (z - 2) + "\t" + nf.format(if (p > 1) 1 else p) + "\t" + (if (p <= pValueThreshold) "present" else "absent"))
          buffer.append(if (p <= pValueThreshold) "1" else "0")
          if (i % 10 == 0)
            buffer.append(" ")

        }
        pw.println("## Binary spoligotype: \n" + buffer.toString())
        val binary = buffer.toString.replaceAll(" ", "").toList.take(43).map(f => f.toString.toInt)
        val octal = binary.grouped(3).map(oc => if (oc.length == 3) (oc(0) * 4 + oc(1) * 2 + oc(2) * 1) else oc(0)).toList
        pw.println("## Octal spoligotype: \n" + octal.mkString(""))
        pw.println("## Lineage:\n" + Lineages.lineage(buffer.toString))
        pw.println()

      }
      pw.close
    } getOrElse {
      println("Could not interpret command-line arguments, quitting!")
      System.exit(-1)
    }

  }
}
