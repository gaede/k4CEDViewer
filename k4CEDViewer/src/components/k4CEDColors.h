#ifndef k4CEDColors_h
#define k4CEDColors_h 1

namespace k4ced{

  /** Class to define the color scheme for the CED event display.
   *  Essentially copied from the CEDViewer package.
   *
   * @author F.Gaede, DESY
   * @date Nov 2025
   */

  class Colors{

  public:
    // access the currently selected colors 
    const std::vector<unsigned>& current() { return _colors ; }

    size_t size() { return _colors.size() ; }
    
    constexpr static int Red         =   0 ;
    constexpr static int Orange      =   1 ;
    constexpr static int Plum        =   2 ;
    constexpr static int Violet      =   3 ;
    constexpr static int Blue        =   4 ;
    constexpr static int LightBlue   =   5 ;
    constexpr static int Aquamarine  =   6 ;
    constexpr static int Green       =   7 ;
    constexpr static int Olive       =   8 ;
    constexpr static int Yellow      =   9 ;
    constexpr static int Dark        =   10 ;
    constexpr static int Light       =   11 ;
    constexpr static int Classic     =   12 ;
    
    constexpr static int ncol = 20 ;
    constexpr static int nscheme = 10 ;

    // use classic if nothing selected
    Colors() : Colors( Classic ) {} 

    // c'tor selects a color scheme
    Colors( unsigned colorScheme ) {

      setColorScheme( colorScheme ) ;
    }


    void setColorScheme( unsigned colScheme ) {

      _colorScheme = colScheme ;

      _colors.clear() ;
      
      switch (_colorScheme) {
	  
        case Red         :
        case Orange      :
        case Plum        :
        case Violet      :
        case Blue        :
        case LightBlue   :
        case Aquamarine  :
        case Green       :
        case Olive       :
        case Yellow      :

            for(int i = 0 ; i<ncol ; ++i)
                _colors.push_back( bojum_colors[  ( _colorScheme * ncol  ) + i ] ) ;
            break;

        case Dark:

            for(int i = 0 ; i<ncol ; ++i )
                for(int j = 0 ; j<nscheme ; ++j )
                    _colors.push_back( bojum_colors[  ( j * ncol  ) + i ] ) ;
            break ;

        case Light:

            for(int i = ncol-1 ; i>=0 ; --i )
                for(int j = 0 ; j<nscheme ; ++j )
                    _colors.push_back( bojum_colors[  ( j * ncol  ) + i ] ) ;
            break ;

        case Classic:
        default:
            _colors.push_back( 0xff00ff );
            _colors.push_back( 0x00ff00 );
            _colors.push_back( 0xffff00 );
            _colors.push_back( 0x0000ff );
            _colors.push_back( 0xff00ff );
            _colors.push_back( 0x00ffff );
            _colors.push_back( 0xffffff );

            _colors.push_back( 0xff88ff );
            _colors.push_back( 0x88ff88 );
            _colors.push_back( 0xffff88 );
            _colors.push_back( 0x8888ff );
            _colors.push_back( 0xff88ff );
            _colors.push_back( 0x88ffff );
            _colors.push_back( 0xffffff );
            break ;
	}
    }
  protected:
    static constexpr unsigned bojum_colors[ncol*nscheme]={
      0x510000,0x660202,0x720202,0x840202,0x960303,0xa80303,0xbc0303,
      0xce0404,0xe00404,0xef0404,0xff0000,0xfc0505,0xf91111,0xf92222,
      0xfc4949,0xf97777,0xf98e8e,0xf9aeae,0xf7b2b2,0xf7cdcd,
      0x512401,0x662d03,0x7a3602,0x934204,0xa54a04,0xb75103,0xc65803,
      0xd35e04,0xe56604,0xff6e00,0xf4710c,0xf4791a,0xf9842a,0xf98f3e,
      0xf99d57,0xf9a768,0xf9b47f,0xf9bf93,0xf9c9a4,0xf9d2b3,
      0x48014c,0x610266,0x700277,0x93039b,0xb103ba,0xc904d3,0xda04e5,
      0xe604f2,0xfa00ff,0xf00ffc,0xec1df7,0xef2ff9,0xeb42f4,0xec58f4,
      0xed6bf4,0xf486f9,0xf59af9,0xf8b5fc,0xf4c3f7,0xf8d6f9,
      0x2d0251,0x3d026d,0x4a0284,0x550399,0x5f03aa,0x6903bc,0x7102cc,
      0x8004e5,0x9800ff,0x8e0ef7,0x9922f9,0xa134f9,0xa845f9,0xb057f9,
      0xbc70f9,0xbf77f9,0xc98ef9,0xd3a4f9,0xddbbf9,0xecd6ff,
      0x00004f,0x020268,0x03037f,0x030399,0x0303b2,0x0404cc,0x0404e0,
      0x0404ef,0x0000ff,0x0c0cf9,0x1b1bf9,0x2a2af9,0x3939f9,0x4d4df9,
      0x6363f9,0x7272f9,0x8484f9,0x9898f9,0xb3b3f9,0xcacaf9,
      0x01434c,0x025c68,0x027382,0x04899b,0x0397aa,0x03abc1,0x04b9d1,
      0x04cbe5,0x00d8ff,0x0cdef9,0x1bddf7,0x2ae1f9,0x3ddff4,0x59e4f7,
      0x6be9f9,0x7cebf9,0x94f0fc,0xa3edf7,0xb6f2f9,0xc7f4f9,
      0x014443,0x01605f,0x027f7d,0x039996,0x03b2af,0x01c6c3,0x04ddda,
      0x04efeb,0x00ffff,0x09f9f5,0x0ef9f5,0x20f9f6,0x32fcf9,0x3ef9f6,
      0x52f9f7,0x6bf9f7,0x7ff9f7,0x95f9f8,0xb1f9f8,0xcaf9f9,
      0x016001,0x027002,0x027f02,0x029102,0x05aa05,0x05bf05,0x06d306,
      0x04e504,0x00ff00,0x09f909,0x18f918,0x2cf92c,0x43f943,0x52f952,
      0x63f963,0x77f977,0x8bf98b,0x9ff99f,0xb3f9b3,0xcaf9ca,
      0x344701,0x4b6603,0x608202,0x739b04,0x83b203,0x96cc04,0xa7e204,
      0xb1ef04,0xb6ff00,0xbaf713,0xbff725,0xc5f73b,0xcbf751,0xd3f968,
      0xd7f97a,0xd8f48b,0xe2f9a2,0xe1f4ad,0xe7f7bb,0xe9f4c8,
      0x565501,0x727002,0x898702,0xa5a303,0xb7b403,0xd1cd04,0xe2df04,
      0xefeb04,0xffff00,0xf9f509,0xf9f618,0xf9f62a,0xf7f43b,0xf9f64a,
      0xf9f759,0xf9f76b,0xf9f77c,0xf9f88b,0xf9f89f,0xfcfbbd
    };
    
    std::vector<unsigned> _colors ;
    unsigned _colorScheme = 0 ;

  };


} // endif

#endif
