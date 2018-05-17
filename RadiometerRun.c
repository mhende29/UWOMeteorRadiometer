/*
 * ADS1256_test.c:
 *	Very simple program to test the serial port. Expects
 *	the port to be looped back to itself
 *
 */
 
/*
             define from bcm2835.h                       define from Board DVK511
                 3.3V | | 5V               ->                 3.3V | | 5V
    RPI_V2_GPIO_P1_03 | | 5V               ->                  SDA | | 5V 
    RPI_V2_GPIO_P1_05 | | GND              ->                  SCL | | GND
       RPI_GPIO_P1_07 | | RPI_GPIO_P1_08   ->                  IO7 | | TX
                  GND | | RPI_GPIO_P1_10   ->                  GND | | RX
       RPI_GPIO_P1_11 | | RPI_GPIO_P1_12   ->                  IO0 | | IO1
    RPI_V2_GPIO_P1_13 | | GND              ->                  IO2 | | GND
       RPI_GPIO_P1_15 | | RPI_GPIO_P1_16   ->                  IO3 | | IO4
                  VCC | | RPI_GPIO_P1_18   ->                  VCC | | IO5
       RPI_GPIO_P1_19 | | GND              ->                 MOSI | | GND
       RPI_GPIO_P1_21 | | RPI_GPIO_P1_22   ->                 MISO | | IO6
       RPI_GPIO_P1_23 | | RPI_GPIO_P1_24   ->                  SCK | | CE0
                  GND | | RPI_GPIO_P1_26   ->                  GND | | CE1

::if your raspberry Pi is version 1 or rev 1 or rev A
RPI_V2_GPIO_P1_03->RPI_GPIO_P1_03
RPI_V2_GPIO_P1_05->RPI_GPIO_P1_05
RPI_V2_GPIO_P1_13->RPI_GPIO_P1_13
::
*/

#include <bcm2835.h>  
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
//CS    -----   SPICS  
//DIN   -----   MOSI
//DOUT  -----   MISO
//SCLK  -----   SCLK
//DRDY  -----   ctl_IO     data starting
//RST   -----   ctl_IO     reset



#define DRDY    RPI_GPIO_P1_11  //P0
#define RST     RPI_GPIO_P1_12  //P1
#define	SPICS	RPI_GPIO_P1_15	//P3

#define CS_1()  bcm2835_gpio_write(SPICS,HIGH)
#define CS_0()  bcm2835_gpio_write(SPICS,LOW)

#define DRDY_IS_LOW()   ((bcm2835_gpio_lev(DRDY)==0))

#define RST_1()     bcm2835_gpio_write(RST,HIGH);
#define RST_0()     bcm2835_gpio_write(RST,LOW);



/* Unsigned integer types  */
#define uint8_t unsigned char
#define uint16_t unsigned short    
#define uint32_t unsigned long     

//Data size is 2^20 which represents approximately 500 seconds worth of data samples
#define DATA_SIZE 1048576


/* Defining Boolean Types  */
typedef enum {FALSE = 0, TRUE = !FALSE} bool;

typedef struct
{
	uint8_t Channel;			/* The current channel*/
	uint8_t ScanMode;	/*Scanning mode,   0  Single-ended input  8 channel£¬ 1 Differential input  4 channel*/
}ADS1256_VAR_T;

typedef struct 
{
    uint32_t header_size;
    uint32_t file_format_version;
    
    char station_code[6];
    char channel[1];
    
    double station_latitude;
    double station_longitude;
    double station_elevation;
    char instrument_string[64];
    
    uint32_t num_samples;
    uint64_t checksum; 
    
    uint32_t unix_s[DATA_SIZE];
    uint32_t unix_us[DATA_SIZE];
    uint32_t intensity[DATA_SIZE];
}rdm;


/*Register definitions Table 23. Register Map --- ADS1256 datasheet Page 30*/
enum
{
	/*Register address, followed by reset the default values */
	REG_STATUS = 0,	 // x1H
	REG_MUX    = 1,  // 01H
	REG_ADCON  = 2,  // 20H
	REG_DRATE  = 3,  // F0H
	REG_IO     = 4,  // E0H
	REG_OFC0   = 5,  // xxH
	REG_OFC1   = 6,  // xxH
	REG_OFC2   = 7,  // xxH
	REG_FSC0   = 8,  // xxH
	REG_FSC1   = 9,  // xxH
	REG_FSC2   = 10, // xxH
};

/* Command definitions TTable 24. Command Definitions --- ADS1256 datasheet Page 34 */
enum
{
	CMD_WAKEUP  = 0x00,	// Completes SYNC and Exits Standby Mode . 0000   0000 (00h)
	CMD_RDATA   = 0x01, // Read Data ............................. 0000   0001 (01h)
	CMD_RDATAC  = 0x03, // Read Data Continuously ................ 0000   0011 (03h)
	CMD_SDATAC  = 0x0F, // Stop Read Data Continuously ........... 0000   1111 (0Fh)
	CMD_RREG    = 0x10, // Read from REG ....................... rrr 0001 rrrr (1xh)
	CMD_WREG    = 0x50, // Write to REG ........................ rrr 0101 rrrr (5xh)
	CMD_SELFCAL = 0xF0, // Offset and Gain Self-Calibration ...... 1111   0000 (F0h)
	CMD_SELFOCAL= 0xF1, // Offset Self-Calibration ............... 1111   0001 (F1h)
	CMD_SELFGCAL= 0xF2, // Gain Self-Calibration ................. 1111   0010 (F2h)
	CMD_SYSOCAL = 0xF3, // System Offset Calibration ............. 1111   0011 (F3h)
	CMD_SYSGCAL = 0xF4, // System Gain Calibration ............... 1111   0100 (F4h)
	CMD_SYNC    = 0xFC, // Synchronize the A/D Conversion ........ 1111   1100 (FCh)
	CMD_STANDBY = 0xFD, // Begin Standby Mode .................... 1111   1101 (FDh)
	CMD_RESET   = 0xFE, // Reset to Power-Up Values .............. 1111   1110 (FEh)
};

rdm  rad_data,config;
ADS1256_VAR_T g_tADS1256;

void  bsp_DelayUS(uint64_t micros);
void ADS1256_StartScan(uint8_t _ucScanMode);
static void ADS1256_Send8Bit(uint8_t _data);
void ADS1256_CfgADC(uint8_t _gain, uint8_t _drate);
static void ADS1256_DelayDATA(void);
static uint8_t ADS1256_Recive8Bit(void);
static void ADS1256_WriteReg(uint8_t _RegID, uint8_t _RegValue);
static uint8_t ADS1256_ReadReg(uint8_t _RegID);
static void ADS1256_WriteCmd(uint8_t _cmd);
uint8_t ADS1256_ReadChipID(void);
static void ADS1256_SetChannal(uint8_t _ch);
static void ADS1256_SetDiffChannal(uint8_t _ch);
static void ADS1256_WaitDRDY(void);
static int32_t ADS1256_ReadData(void);

uint8_t samplerate(double ratenum);
void Init_ADC(double _gain,double _sps,uint8_t _mode);
int32_t Read_Single_Channel(uint8_t channel);
void Init_Single_Channel(uint8_t channel);
uint8_t getgain(double givengain);
void savedat(rdm *data,const char *path);
int Runtime(double time);


/*
*********************************************************************************************************
*	Name: bsp_DelayUS
*	Description: A delay of a desired amount of microseconds
*	Arguments: 
*       micros: The desired amount of microseconds
*	Return: NULL
*********************************************************************************************************
*/

void  bsp_DelayUS(uint64_t micros)
{
		bcm2835_delayMicroseconds (micros);
}

/*
*********************************************************************************************************
*	Name: bsp_InitADS1256
*	Description: Configuration of the STM32 GPIO and SPI interface£¬The connection ADS1256
*	Arguments: NULL
*	Return: NULL
*********************************************************************************************************
*/

void bsp_InitADS1256(void)
{
#ifdef SOFT_SPI
	CS_1();
	SCK_0();
	DI_0();
#endif

}

/*
*********************************************************************************************************
*	Name: ADS1256_StartScan
*	Description: Sets the scan mode to single input or differential input
*	Arguments: 
*       _ucScanMode : 0 Single-ended input, 1 Differential input
*	Return: NULL
*********************************************************************************************************
*/

void ADS1256_StartScan(uint8_t _ucScanMode)
{
	g_tADS1256.ScanMode = _ucScanMode;
	
	{
		uint8_t i;

		g_tADS1256.Channel = 0;
	}

}

/*
*********************************************************************************************************
*	Name: ADS1256_Send8Bit
*	Description: SPI bus to send 8 bit data
*	Arguments: 
*       _data:  The data to be sent
*	Return: NULL
*********************************************************************************************************
*/

static void ADS1256_Send8Bit(uint8_t _data)
{

	bsp_DelayUS(2);
	bcm2835_spi_transfer(_data);
}

/*
*********************************************************************************************************
*	Name: ADS1256_CfgADC
*	Description: The configuration parameters of ADC, gain and data rate
*	Arguments: 
*       _gain:  Desired gain 
*       _drate: Desired samples per second
*	Return: NULL
*********************************************************************************************************
*/

void ADS1256_CfgADC(uint8_t _gain, uint8_t _drate)
{

    // Wait before sending data
	ADS1256_WaitDRDY();

	{
        /* Storage ads1256 register configuration parameters */
        uint8_t buf[4];		

		/*Status register define
			Bits 7-4 ID3, ID2, ID1, ID0  Factory Programmed Identification Bits (Read Only)

			Bit 3 ORDER: Data Output Bit Order
				0 = Most Significant Bit First (default)
				1 = Least Significant Bit First
			Input data  is always shifted in most significant byte and bit first. Output data is always shifted out most significant
			byte first. The ORDER bit only controls the bit order of the output data within the byte.

			Bit 2 ACAL : Auto-Calibration
				0 = Auto-Calibration Disabled (default)
				1 = Auto-Calibration Enabled
			When Auto-Calibration is enabled, self-calibration begins at the completion of the WREG command that changes
			the PGA (bits 0-2 of ADCON register), DR (bits 7-0 in the DRATE register) or BUFEN (bit 1 in the STATUS register)
			values.

			Bit 1 BUFEN: Analog Input Buffer Enable
				0 = Buffer Disabled (default)
				1 = Buffer Enabled

			Bit 0 DRDY :  Data Ready (Read Only)
				This bit duplicates the state of the DRDY pin.

			ACAL=1  enable  calibration
		*/
		//buf[0] = (0 << 3) | (1 << 2) | (1 << 1);//enable the internal buffer
        buf[0] = (0 << 3) | (1 << 2) | (0 << 1);  // The internal buffer is prohibited

        //ADS1256_WriteReg(REG_STATUS, (0 << 3) | (1 << 2) | (1 << 1));

		buf[1] = 0x08;	

		/*	ADCON: A/D Control Register (Address 02h)
			Bit 7 Reserved, always 0 (Read Only)
			Bits 6-5 CLK1, CLK0 : D0/CLKOUT Clock Out Rate Setting
				00 = Clock Out OFF
				01 = Clock Out Frequency = fCLKIN (default)
				10 = Clock Out Frequency = fCLKIN/2
				11 = Clock Out Frequency = fCLKIN/4
				When not using CLKOUT, it is recommended that it be turned off. These bits can only be reset using the RESET pin.

			Bits 4-3 SDCS1, SCDS0: Sensor Detect Current Sources
				00 = Sensor Detect OFF (default)
				01 = Sensor Detect Current = 0.5 ŠÌ A
				10 = Sensor Detect Current = 2 ŠÌ A
				11 = Sensor Detect Current = 10ŠÌ A
				The Sensor Detect Current Sources can be activated to verify  the integrity of an external sensor supplying a signal to the
				ADS1255/6. A shorted sensor produces a very small signal while an open-circuit sensor produces a very large signal.

			Bits 2-0 PGA2, PGA1, PGA0: Programmable Gain Amplifier Setting
				000 = 1 (default)
				001 = 2
				010 = 4
				011 = 8
				100 = 16
				101 = 32
				110 = 64
				111 = 64
		*/
		buf[2] = (0 << 5) | (0 << 3) | (_gain << 0);
		//ADS1256_WriteReg(REG_ADCON, (0 << 5) | (0 << 2) | (GAIN_1 << 1));	/*choose 1: gain 1 ;input 5V/
		buf[3] = _drate;	

        // Send the configuration parameters to the ADC
		CS_0();	
		ADS1256_Send8Bit(CMD_WREG | 0);	/* Write command register, send the register address */
		ADS1256_Send8Bit(0x03);			/* Register number 4,Initialize the number  -1*/

		ADS1256_Send8Bit(buf[0]);	/* Set the status register */
		ADS1256_Send8Bit(buf[1]);	/* Set the input channel parameters */
		ADS1256_Send8Bit(buf[2]);	/* Set the ADCON control register,gain */
		ADS1256_Send8Bit(buf[3]);	/* Set the output rate */

		CS_1();	/* SPI  cs = 1 */
	}

	bsp_DelayUS(50);
}

/*
*********************************************************************************************************
*	Name: ADS1256_DelayDATA
*	Description: Delay of 10 microseconds
*	Arguments: NULL
*	Return: NULL
*********************************************************************************************************
*/

static void ADS1256_DelayDATA(void)
{
	/*
		Delay from last SCLK edge for DIN to first SCLK rising edge for DOUT: RDATA, RDATAC,RREG Commands
		min  50   CLK = 50 * 0.13uS = 6.5uS
	*/
	bsp_DelayUS(10);	/* The minimum time delay 6.5us */
}

/*
*********************************************************************************************************
*	Name: ADS1256_Recive8Bit
*	Description: SPI bus receive function
*	Arguments: NULL
*	Return: 
*       read: The byte that is read in from the bus
*********************************************************************************************************
*/

static uint8_t ADS1256_Recive8Bit(void)
{
	uint8_t read = 0;
	read = bcm2835_spi_transfer(0xff);
	return read;
}

/*
*********************************************************************************************************
*	Name: ADS1256_WriteReg
*	Description: Write the corresponding register
*	Arguments: 
*       _RegID: Register ID
*		_RegValue: Register Value
*	Return: NULL
*********************************************************************************************************
*/

static void ADS1256_WriteReg(uint8_t _RegID, uint8_t _RegValue)
{
	CS_0();	                                /*        SPI  cs = 0        */
	ADS1256_Send8Bit(CMD_WREG | _RegID);	/*   Write command register  */
	ADS1256_Send8Bit(0x00);		            /* Write the register number */
                                            /*                           */
	ADS1256_Send8Bit(_RegValue);	        /*    Send register value    */
	CS_1();	                                /*        SPI  cs = 1        */
}

/*
*********************************************************************************************************
*	Name: ADS1256_ReadReg
*	Description: Read the corresponding register
*	Arguments: 
*       _RegID: Register ID
*	Return: 
*       read: Register value read
*********************************************************************************************************
*/

static uint8_t ADS1256_ReadReg(uint8_t _RegID)
{
	uint8_t read;

	CS_0();	                                /*         SPI cs = 0        */
	ADS1256_Send8Bit(CMD_RREG | _RegID);	/*  Write command register   */
	ADS1256_Send8Bit(0x00);	                /* Write the register number */
                                            /*                           */
	ADS1256_DelayDATA();	                /*         delay time        */
                                            /*                           */
	read = ADS1256_Recive8Bit();	        /*  Read the register values */
	CS_1();	                                /*         SPI cs = 1        */

	return read;
}

/*
*********************************************************************************************************
*	Name: ADS1256_WriteCmd
*	Description: Sending a single byte order
*	Arguments: 
*       _cmd : Desired command
*	Return: NULL
*********************************************************************************************************
*/

static void ADS1256_WriteCmd(uint8_t _cmd)
{
	CS_0();	                /* SPI   cs = 0 */
	ADS1256_Send8Bit(_cmd); /*              */
	CS_1();	                /* SPI  cs  = 1 */
}

/*
*********************************************************************************************************
*	Name: ADS1256_ReadChipID
*	Description: Read the chip ID
*	Arguments: NULL
*	Return: 
*       id: status register shift 4 bits
*********************************************************************************************************
*/

uint8_t ADS1256_ReadChipID(void)
{
	uint8_t id;

	ADS1256_WaitDRDY();
	id = ADS1256_ReadReg(REG_STATUS);
	return (id >> 4);
}

/*
*********************************************************************************************************
*	Name: ADS1256_SetChannal
*	Description: Configuration channel number
*	Arguments:  
*       _ch: Channel number 0--7
*	Return: NULL
*********************************************************************************************************
*/

static void ADS1256_SetChannal(uint8_t _ch)
{
	/*
	Bits 7-4 PSEL3, PSEL2, PSEL1, PSEL0: Positive Input Channel (AINP) Select
		0000 = AIN0 (default)
		0001 = AIN1
		0010 = AIN2 (ADS1256 only)
		0011 = AIN3 (ADS1256 only)
		0100 = AIN4 (ADS1256 only)
		0101 = AIN5 (ADS1256 only)
		0110 = AIN6 (ADS1256 only)
		0111 = AIN7 (ADS1256 only)
		1xxx = AINCOM (when PSEL3 = 1, PSEL2, PSEL1, PSEL0 are ¡°don¡¯t care¡±)

		NOTE: When using an ADS1255 make sure to only select the available inputs.

	Bits 3-0 NSEL3, NSEL2, NSEL1, NSEL0: Negative Input Channel (AINN)Select
		0000 = AIN0
		0001 = AIN1 (default)
		0010 = AIN2 (ADS1256 only)
		0011 = AIN3 (ADS1256 only)
		0100 = AIN4 (ADS1256 only)
		0101 = AIN5 (ADS1256 only)
		0110 = AIN6 (ADS1256 only)
		0111 = AIN7 (ADS1256 only)
		1xxx = AINCOM (when NSEL3 = 1, NSEL2, NSEL1, NSEL0 are ¡°don¡¯t care¡±)
	*/
	if (_ch > 7)
	{
		return;
	}
    ADS1256_WriteReg(REG_MUX,0x01);
	ADS1256_WriteReg(REG_MUX, (_ch << 4) | (1 << 3));	/* Bit3 = 1, AINN connection AINCOM */
}

/*
*********************************************************************************************************
*	Name: ADS1256_SetDiffChannal
*	Description: The configuration difference channel
*	Arguments:  
*       _ch:  Channel number 0--3
*	Return: NULL
*********************************************************************************************************
*/

static void ADS1256_SetDiffChannal(uint8_t _ch)
{
	/*
	Bits 7-4 PSEL3, PSEL2, PSEL1, PSEL0: Positive Input Channel (AINP) Select
		0000 = AIN0 (default)
		0001 = AIN1
		0010 = AIN2 (ADS1256 only)
		0011 = AIN3 (ADS1256 only)
		0100 = AIN4 (ADS1256 only)
		0101 = AIN5 (ADS1256 only)
		0110 = AIN6 (ADS1256 only)
		0111 = AIN7 (ADS1256 only)
		1xxx = AINCOM (when PSEL3 = 1, PSEL2, PSEL1, PSEL0 are ¡°don¡¯t care¡±)

		NOTE: When using an ADS1255 make sure to only select the available inputs.

	Bits 3-0 NSEL3, NSEL2, NSEL1, NSEL0: Negative Input Channel (AINN)Select
		0000 = AIN0
		0001 = AIN1 (default)
		0010 = AIN2 (ADS1256 only)
		0011 = AIN3 (ADS1256 only)
		0100 = AIN4 (ADS1256 only)
		0101 = AIN5 (ADS1256 only)
		0110 = AIN6 (ADS1256 only)
		0111 = AIN7 (ADS1256 only)
		1xxx = AINCOM (when NSEL3 = 1, NSEL2, NSEL1, NSEL0 are ¡°don¡¯t care¡±)
	*/
	if (_ch == 0)
	{
		ADS1256_WriteReg(REG_MUX, (0 << 4) | 1);	/* DiffChannal  AIN0-AIN1 */
	}
	else if (_ch == 1)
	{
		ADS1256_WriteReg(REG_MUX, (2 << 4) | 3);	/*DiffChannal   AIN2-AIN3 */
	}
	else if (_ch == 2)
	{
		ADS1256_WriteReg(REG_MUX, (4 << 4) | 5);	/*DiffChannal    AIN4--AIN5 */
	}
	else if (_ch == 3)
	{
		ADS1256_WriteReg(REG_MUX, (6 << 4) | 7);	/*DiffChannal   AIN6-AIN7 */
	}
}

/*
*********************************************************************************************************
*	Name: ADS1256_WaitDRDY
*	Description: Delay time  wait for automatic calibration
*	Arguments: NULL
*	Return: NULL
*********************************************************************************************************
*/

static void ADS1256_WaitDRDY(void)
{
	uint32_t i;

	for (i = 0; i < 400000; i++)
	{
		if (DRDY_IS_LOW())
		{
			break;
		}
	}
	if (i >= 400000)
	{
		printf("ADS1256_WaitDRDY() Time Out ...\r\n");		
	}
}

/*
*********************************************************************************************************
*	Name: ADS1256_ReadData
*	Description: Read ADC value
*	Arguments: NULL
*	Return: 
*      read: Return the 24 bit number read by the ADC
*********************************************************************************************************
*/

static int32_t ADS1256_ReadData(void)
{
	uint32_t read = 0;
    static uint8_t buf[3];

	CS_0();	
    
    /* read ADC command */
	ADS1256_Send8Bit(CMD_RDATA);	

    /* delay time */
	ADS1256_DelayDATA();	

	/* Read the sample results 24bit 1 byte at a time */
    buf[0] = ADS1256_Recive8Bit();
    buf[1] = ADS1256_Recive8Bit();
    buf[2] = ADS1256_Recive8Bit();

    /* Combine the 3 read bytes into one 24 bit number */
    read = ((uint32_t)buf[0] << 16) & 0x00FF0000;
    read |= ((uint32_t)buf[1] << 8);  
    read |= buf[2];

	CS_1();	

	/* Extend a signed number*/
    if (read & 0x800000)
    {
	    read |= 0xFF000000;
    }

	return (int32_t)read;
}

/*
*********************************************************************************************************
*	Name: Init_ADC
*	Description:  Initializes the micro and sets the gain, sample rate and mode of the ADC
*	Arguments: 
*       _gain: Is the gain
*       _sps : Is the sample rate
*       _mode: Is the mode, either single or differential
*	Return: NULL
*********************************************************************************************************
*/

void Init_ADC(double _gain,double _sps,uint8_t _mode)
{
	
    bcm2835_init();
    bcm2835_spi_begin();
    bcm2835_spi_setBitOrder(BCM2835_SPI_BIT_ORDER_LSBFIRST );    // The default
    bcm2835_spi_setDataMode(BCM2835_SPI_MODE1);                  // The default
    bcm2835_spi_setClockDivider(BCM2835_SPI_CLOCK_DIVIDER_1024); // The default
    bcm2835_gpio_fsel(SPICS, BCM2835_GPIO_FSEL_OUTP);
    bcm2835_gpio_write(SPICS, HIGH);
    bcm2835_gpio_fsel(DRDY, BCM2835_GPIO_FSEL_INPT);
    bcm2835_gpio_set_pud(DRDY, BCM2835_GPIO_PUD_UP);
    
    uint8_t gain=getgain(_gain);
    uint8_t sps=samplerate(_sps);
    
    ADS1256_CfgADC(gain, sps);
    ADS1256_StartScan(_mode);
}

/*
*********************************************************************************************************
*	Name: samplerate
*	Description: Gives the hex value associated to the sample rate
*	Arguments: 
*       ratenum : The decimal version of the rate	   
*	Return:  
*       rate : The samples per second hex code 
*********************************************************************************************************
*/

uint8_t samplerate(double ratenum)
{
    uint8_t rate;
    if(ratenum==30000){
        rate=0xF0;
    }
    else if (ratenum==15000){
        rate=0xE0;
    }
    else if (ratenum==7500){
        rate=0xD0;
    }
    else if (ratenum==3750){
        rate=0xC0;
    }
    else if (ratenum==2000){
        rate=0xB0;
    }
    else if (ratenum==1000){
        rate=0xA1;
    }
    else if (ratenum==500){
        rate=0x92;
    }
    else if (ratenum==100){
        rate=0x82;
    }
    else if (ratenum==60){
        rate=0x72;
    }
    else if (ratenum==50){
        rate=0x63;
    }
    else if (ratenum==30){
        rate=0x53;
    }
    else if (ratenum==25){
        rate=0x43;
    }
    else if (ratenum==15){
        rate=0x33;
    }
    else if (ratenum==10){
        rate=0x23;
    }
    else if (ratenum=5){
        rate=0x13;
    }
    else if (ratenum==2.5){
        rate=0x03;
    }
    else{
        rate=0xF0;
    }
    return rate;
} 

/*
*********************************************************************************************************
*	Name: getgain
*	Description: Gives the value associated to the given gain
*	Arguments: 
*       givengain : The decimal value of the gain
*	Return:  
*       gain : The gain code 
*********************************************************************************************************
*/

uint8_t getgain(double givengain)
{
    uint8_t gain;
    
    if(givengain==1){
        gain=0;
    }
    else if (givengain==2){
        gain=1;
    }
    else if (givengain==4){
        gain=2;
    }
    else if (givengain==8){
        gain=3;
    }
    else if (givengain==16){
        gain=4;
    }
    else if (givengain==32){
        gain=5;
    }
    else if (givengain==64){
        gain=6;
    }
    else{
        gain=0;
    }
    
    return gain;
}

/*
*********************************************************************************************************
*	Name: Read_Single_Channel
*	Description: Reads a given channel
*	Arguments: 
*       channel: The channel number
*	Return: 
*       The ADU value read
*********************************************************************************************************
*/

int32_t Read_Single_Channel(uint8_t channel)
{
    // 0 means not done, 1 means done
    uint8_t done=0; 
    
    // Confirms the channel we want to read is the one were set to
    if(g_tADS1256.Channel!=channel)
    {
        // If it's not the desired channel, switch to it
        Init_Single_Channel(channel);
    }
    
    while(done == 0){
        
        // While the data ready pin is waiting for a new sample loop 
        // until its ready
        if (DRDY_IS_LOW())
        {
            // When the data is ready sync up
            ADS1256_WriteCmd(CMD_SYNC);
            bsp_DelayUS(5);
    
            // Complete the sync
            ADS1256_WriteCmd(CMD_WAKEUP);
            bsp_DelayUS(25);
            
            // Acknowledge data is ready to be read
            done=1;
        }
    }
    
    // Return the sampled value
    return ADS1256_ReadData();	
}

/*
*********************************************************************************************************
*	Name: savedat
*	Description: Saves the rdm struct to the desired file path
*	Arguments: 
*       *data : Input the address of the struct. Ex: &my_struct
*       *path : Input the desired path string. Ex: "//home//path_to_desired_file//"
*        Note : Must use two // in order to properly format the string
*	Return: NULL
*********************************************************************************************************
*/

void savedat(rdm *data, const char *path)
{
    // Define length of filename including the path
    char buf[128];                                   
    
    // CCNNNN_X_YYYYMMDD-HHMMSS.UUUUUU_hhmmss.uuuuuu.rdm is the full extension CCNNNN is the station code string
    // Get the time_t version of the first sample in the data structure
    time_t t1 = (time_t)data->unix_s[0];           
    
    // Get the time_t version of the final sample in the data structure   
    time_t t2 = (time_t)data->unix_s[DATA_SIZE-1];    
    
    // jcal = julian calendar variables
    struct tm jcal1, jcal2;           
    
    // Get the julian equivalent of the unix time               
    gmtime_r(&t1, &jcal1);  
    
    // Get the julian equivalent of the unix time                          
    gmtime_r(&t2, &jcal2);                            
                                                   
    // The conversion of all this data into a string is done below      
    // %.*s prints out a defined number of characters from a provided string
    snprintf(buf,127,"%s%.*s_%s_%04d%02d%02d-%02d%02d%02d.%06d_%02d%02d%02d.%06d.rdm", path, 6,
    data->station_code, data->channel, jcal1.tm_year+1900, jcal1.tm_mon+1, jcal1.tm_mday, jcal1.tm_hour,
    jcal1.tm_min, jcal1.tm_sec, data->unix_us[0], jcal2.tm_hour, jcal2.tm_min, jcal2.tm_sec,
    data->unix_us[DATA_SIZE-1]);
    
    // Open a file and save the data in binary
    FILE *file_path;
    file_path=fopen(buf, "wb");
    fwrite(data,sizeof(*data), 1, file_path);
    
}

/*
*********************************************************************************************************
*	Name: Init_Single_Channel
*	Description: Initializes a channel
*	Arguments: The channel
*	Return:  NULL
*********************************************************************************************************
*/
void Init_Single_Channel(uint8_t channel)
{
    // Record the current channel being sampled
    g_tADS1256.Channel = channel;
    ADS1256_SetChannal(channel);	/*Switch channel mode */
    bsp_DelayUS(5);
}

/*
*********************************************************************************************************
*	Name: Init_Single_Channel
*	Description: Initializes a channel
*	Arguments: The channel
*	Return:  NULL
*********************************************************************************************************
*/
int Runtime(double time)
{
    int run;
    if(time==-1){
        run=-1;
        return run;
    }
    else{
        // One run is 500s and one hour has 3600s, therefore the ratio of 36/5 is how many runs can be done 
        // in one hour. 36/5=7.2
        run=(int)ceil(time*(7.2));
        return run;
    }
}

/*
*********************************************************************************************************
*	Name: main
*	Description:  
*	Arguments: NULL
*	Return:  NULL
*********************************************************************************************************
*/

int  main()
{
    uint32_t num =0;
  	int32_t adc=0;
    uint64_t cksum=0;
    uint32_t count=0;
    
    int loops =0;
    int i=0;
    
    struct timespec tp; 
    
    
    //// Configuration /////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Number of hours to run, -1 means run forever
    double run_hrs = -1;
    
    //Mode is either 0 for single input or 1 for differential
    uint8_t mode = 0;
    
    //Gains are 1,2,4,8,16,32 and 64
    double gain = 1; 
    
    // Function sample rate determines the equivalent hex value associated to the sps
    // Inputs must be 30000,15000,7500,3750,2000,1000,500,100,60,50,30,25,15,10,5,2.5
    double sps = 3750; 
    
    // 0 for one channel, 1 for multiple channels
    uint8_t single_channel=0; 

    uint8_t channel=0;
    
    // Define a configuration for the radiometer
    config=(rdm){
        
        // Size of header in bytes including 
        .header_size=120,
        
        
        .file_format_version=1,
        
        
        .station_code="CA0001",
        
        
        .channel="A",
        
        
        .station_latitude=43.19279,
        
        
        .station_longitude=-81.31566,
        
        
        .station_elevation=324.0,
        
        
        .instrument_string="Radiometer prototype.",
        
        };
        
        
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    loops = Runtime(run_hrs);
    
    Init_ADC(gain,sps,mode);
    Init_Single_Channel(channel);
	
    
    while(1){
        
        printf("Recording a new file");
        rad_data=config;
        
        while(count<DATA_SIZE){
        
            adc= Read_Single_Channel(channel);
        
            clock_gettime(CLOCK_REALTIME,&tp);
        
            cksum+=adc;
        
            rad_data.intensity[count]=adc;
        
            rad_data.unix_s[count]=tp.tv_sec;
        
            rad_data.unix_us[count]=tp.tv_nsec/1000;
        
            count++;
        }   
            
        rad_data.num_samples=count;
        rad_data.checksum=cksum;
        savedat(&rad_data,"//home//pi//Desktop//Michael//RDM//");
        rad_data=(rdm){0};
        count=0;
        if(loops!=-1){
            i++;
            if(i==loops){
                break;
            }
        }
    }

    bcm2835_spi_end();
    bcm2835_close();
    return 0;
}

// printf("%lld\n",(long long)rad_data.checksum); // How to print long numbers
//if(count==1000){
//                count=DATA_SIZE-1;
//            }
