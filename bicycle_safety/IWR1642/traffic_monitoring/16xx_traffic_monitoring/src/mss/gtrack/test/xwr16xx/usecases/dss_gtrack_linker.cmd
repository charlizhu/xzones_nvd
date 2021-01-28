/*----------------------------------------------------------------------------*/
/* Linker Settings                                                            */
--retain="*(.intvecs)"

-stack 0x1000

/*----------------------------------------------------------------------------*/
/* Section Configuration                                                      */
SECTIONS
{
/*    systemHeap : {} > 0x007E1000*/
    systemHeap : {} > 0x00800000	
}
/*----------------------------------------------------------------------------*/

