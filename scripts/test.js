var fs = require('fs');
var bam_input_file = '/mnt/e/2018_Xuyu_project/BAMS/BRCA_RNAseq_Realign/1f578ca2-3646-44f7-a1d3-76491344ac27_gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
var readStream = fs.createReadStream(bam_input_file);
var a = 0;
readStream
.on('readable', () => {
  if (a == 0) {
    a = 1;
    var chunk;
    while (null !== (chunk = readStream.read())) {
      console.log(chunk.length);
      console.log('  ID1:   ' + chunk.readUInt8(0)); // 1 byte
      console.log('  ID2:   ' + chunk.readUInt8(1)); // 2 bytes
      console.log('  CM:    ' + chunk.readUInt8(2)); // 3 bytes
      console.log('  FLG:   ' + chunk.readUInt8(3)); // 4 bytes
      console.log('  MTIME: ' + chunk.readUInt32LE(4));
      console.log('  XFL:   ' + chunk.readUInt8(8));
      console.log('  OS:    ' + chunk.readUInt8(9));
      console.log('  XLEN:  ' + chunk.readUInt16LE(10)); // NUMBER OF BYTES in extra
      console.log('    SI1:   ' + chunk.readUInt8(12));
      console.log('    SI2:   ' + chunk.readUInt8(13));
      console.log('    SLEN:  ' + chunk.readUInt16LE(14));
      console.log('    BSIZE: ' + chunk.readUInt16LE(16));
      // console.log('  CDATA: ' + chunk.readUInt8(18));
      // console.log('  CDATA: ' + chunk.toString('ascii',18, 20));
      console.log('  CRC32: ' + chunk.readUInt32LE(702));
      console.log('  ISIZE: ' + chunk.readUInt32LE(706));
      console.log('READS');
      for (var i = 16; i < 500; ++i) {
        console.log('  CDATA: ' + chunk.toString('ascii',i, i+4));
      }
      console.log(' magic: ' + chunk.toString('ASCII', 710, 715));
    }
  }
})
.on('end', () => {
  console.log('ENDED');
});
