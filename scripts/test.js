var fs = require('fs');
var zlib = require('zlib');
var refseq = '=ACMGRSVTWYHKDBN';
var cigmap = 'MIDNSHP=X';

function processBGZFHeader(chunk) {
  console.log('-------------------------------')
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
  var block = chunk.readUInt16LE(16) - chunk.readUInt16LE(10)
  // console.log('  CDATA: ' + chunk.readUInt8(18));
  var buffer = chunk.slice(18, block);
  zlib.inflateRaw(buffer, function(err, buffer) {
    if (!err) {
      console.log('-------------------------------')
      console.log(buffer.length);
      console.log('    magic:   ' + buffer.toString('ascii', 0, 5)); // 1 byte
      console.log('    l_text:   ' + buffer.readInt32LE(4)); // 1 byte
      var l_text = buffer.readInt32LE(4);
      console.log('    header:\n   ' + buffer.toString('ascii', 8, l_text + 8)); // 1 byte
      console.log('    n_ref:   ' + buffer.readInt32LE(l_text + 8)); // 1 byte
      console.log('    l_name:   ' + buffer.readInt32LE(l_text + 8 + 4)); // 1 byte
      var l_name = buffer.readInt32LE(l_text + 8 + 4);
      console.log('    name:   ' + buffer.toString('ascii', l_text + 16, l_name + l_text + 16)); // 1 byte
      console.log('    l_ref:   ' + buffer.readInt32LE(l_name + l_text + 16)); // 1 byte
      // console.log('    b_size:   ' + buffer.readInt32LE(2421)); // 1 byte
      console.log('-------------------------------')
    }else{
      console.log('\ninflate error:\n'+err);
    }
  });
  console.log(block);
  console.log('  CRC32: ' + chunk.readUInt32LE(block-1));
  console.log('  ISIZE: ' + chunk.readUInt32LE(block+3));
  console.log('-------------------------------')
  return chunk.readUInt16LE(16) + 1;
}


function processBGZFBody(chunk) {
  console.log('-------------------------------')
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
  var block = chunk.readUInt16LE(16) - chunk.readUInt16LE(10)
  // console.log('  CDATA: ' + chunk.readUInt8(18));
  var buffer = chunk.slice(18, block);
  zlib.inflateRaw(buffer, function(err, buffer) {
    if (!err) {
      console.log('-------------------------------')
      console.log(buffer.length);
      for (var i = 0; i < 2; ++i) {
        console.log('    b_size:   ' + buffer.readInt32LE(0)); // 1 byte
        console.log('    refid:   ' + buffer.readInt32LE(4)); // 1 byte
        console.log('    pos:   ' + buffer.readInt32LE(8)); // 1 byte
        console.log('    lreadname:   ' + buffer.readUInt8(12)); // 1 byte
        console.log('    mapq:   ' + buffer.readUInt8(13)); // 1 byte
        console.log('    bin:   ' + buffer.readUInt16LE(14)); // 1 byte
        console.log('    n_cigar_op:   ' + buffer.readUInt16LE(16)); // 1 byte
        console.log('    flag:   ' + buffer.readUInt16LE(18)); // 1 byte
        console.log('    l_seq:   ' + buffer.readInt32LE(20)); // 1 byte
        console.log('    nrefid:   ' + buffer.readInt32LE(24)); // 1 byte
        console.log('    npos:   ' + buffer.readInt32LE(28)); // 1 byte
        console.log('    tlen:   ' + buffer.readInt32LE(32)); // 1 byte
        // TODO: add check for null character
        console.log('    rname:   ' + buffer.toString('ascii', 36, 36 + buffer.readUInt8(12))); // 1 byte
        var pos = 36 + buffer.readUInt8(12);
        for (var cig = 0; cig < buffer.readUInt16LE(16); cig++) {

          //012345678
          var tst =  buffer.readUInt32LE(pos);
          var op = tst & 0xF;
          var op_len = tst >>> 4;
          console.log('cig : ' + cigmap[op] + '' + op_len);
          pos += 4;
        }
        var char_seq = ''
        for (var char = 0; char < Math.floor((buffer.readInt32LE(20)+1)/2); char++) {
          var num = buffer.readUInt8(pos);
          var a = num >> 4 & 0xF;
          var b = num & 0xF;
          char_seq += refseq[a];
          char_seq += refseq[b];
          // console.log('    seq:   ' + num + '-' + a + ' : ' + b); // 1 byte
          pos += 1;
        }
        console.log('    seq:    ' + char_seq);
        console.log('   seq_len: ' + char_seq.length);
        console.log(pos);
        qual = ''
        for (var char = 0; char < buffer.readInt32LE(20); char++) {
            qual += String.fromCharCode(buffer.readUInt8(pos + char) + 33);
        }
        console.log('    qual:   |' + qual + '|');// + buffer.readInt32LE(20)+'|')); // 1 byte
        pos += buffer.readInt32LE(20)+1;
        while (pos < buffer.readInt32LE(0) + 2) {
          console.log('    tag:   |' + buffer.toString('ascii', pos, pos + 2)+'|'); // 1 byte
          pos += 2;
          console.log('    val_type:   ' + buffer.toString('ascii', pos, pos + 1)); // 1 byte
          var val_type = buffer.toString('ascii', pos, pos + 1);
          pos += 1;
          if (val_type == 'C') {
            console.log('    val:   ' + buffer.readUInt8(pos)); // 1 byte
            pos += 1;
          } else if (val_type == 'Z'){
            var outp = ''
            while (true) {

              if (buffer.toString('ascii', pos, pos+1) == '\0') {
                break;
              }
              outp += buffer.toString('ascii', pos, pos+1);
              pos += 1
            }
            pos += 1;
            console.log('    val:   |' + outp + '|'); // 1 byte
          } else {
            // AcCsSiIfZHB
            console.log('uncoded val type: ' + val_type + ' pos ' + pos);
            pos += 1;
          }
        }
        console.log('---------------------' + pos);
        buffer = buffer.slice(pos, buffer.length);
      }
    }else{
      console.log('\ninflate error:\n'+err);
    }
  });
  console.log(block);
  console.log('  CRC32: ' + chunk.readUInt32LE(block-1));
  console.log('  ISIZE: ' + chunk.readUInt32LE(block+3));
  console.log('-------------------------------')
  return chunk.readUInt16LE(16) + 1;
}


var bam_input_file = '../test_data/3e25dd86-256f-4b4a-bd54-8d8e83d47e37_gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
var bam_index_file = '../test_data/3e25dd86-256f-4b4a-bd54-8d8e83d47e37_gdc_realn_rehead.Aligned.sortedByCoord.out.bam.bai'
var readStream = fs.createReadStream(bam_input_file);
var a = 0;
readStream
.on('readable', () => {
  if (a == 0) {
    a = 1;
    var chunk;
    while (null !== (chunk = readStream.read())) {

      //These are for BAM
      console.log(chunk.length);
      var offset = 0
      for (var i = 0; i < 1;++i) {

          offset = processBGZFHeader(chunk.slice(offset, chunk.length));
          console.log(offset);
      }
      offset = processBGZFBody(chunk.slice(offset, chunk.length));
      console.log(offset);
    }
  }
})
.on('end', () => {
  console.log('ENDED');
});
