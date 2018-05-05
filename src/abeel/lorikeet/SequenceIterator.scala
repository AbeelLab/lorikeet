package abeel.lorikeet

import java.io.{File, FileInputStream}

import atk.compbio.fastq.FastQFile
import atk.compbio.fastq.FastQFile.FastqRecord
import net.sf.samtools.SAMFileReader
import java.util.zip.GZIPInputStream

object SequenceIterator {
  class SAMIterator(file: File) extends Iterator[FastqRecord]  {
    private val reader = new SAMFileReader(file)
    private val it = reader.iterator()

    def hasNext(): Boolean = it.hasNext

    def next(): FastqRecord = {
      val sr = it.next()
      FastqRecord(sr.getReadName, sr.getReadString, sr.getBaseQualityString)
    }
  }

  object GzippedFastqFile {
    def apply(file: File): Iterator[FastqRecord] = {
      val in = new GZIPInputStream(new FileInputStream(file))
      io.Source.fromInputStream(in).getLines.filterNot(_.trim.length == 0).grouped(4).map(
        group => { FastqRecord(group(0).drop(1), group(1), group(3)) }
      )
    }
  }

  def apply(file: File): Iterator[FastqRecord] = {
    if (file.getPath.endsWith(".sam") || file.getPath.endsWith(".bam"))
      new SAMIterator(file)
    else if (file.getPath.endsWith(".gz"))
      GzippedFastqFile(file)
    else
      FastQFile(file)
  }
}