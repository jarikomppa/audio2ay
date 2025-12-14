        DEVICE ZXSPECTRUM48         ; Device setting for sjasmplus (.tap writing etc)

        ORG $8000                   ; Let's start our code at 32k
        di                          ; Disable interrupts
        ld  sp,     0x8000          ; Set stack to grow down from our code
        ei

        ; Set up (otherwise untouched) audio regs

        ld bc, 0xFFFD 
        ld a, 7 ; AY mixer
        out (c), a
        ld bc, 0xBFFD
        ld a, 8+16+32; 111000 for tone on, noise off on all channels
        out (c), a

        ld hl, data+4 ; start of data
        ld de, (data+2) ; number of blocks
        ld (datacount), de
        
mainloop:
        ld a, (data) ; number of AY channels (1-3)
        ld b, a
ayloop:
        push bc
        ld e, b
        ld a, b ; 1,2,3
        dec a ; 0,1,2
        add a,a ; 0,2,4
        ld bc, 0xFFFD 
        out (c), a ; fine tune for channel (0,2,4)
        ld d, a
        ld a, (hl)
        ld bc, 0xBFFD
        out (c), a
        inc hl
        ld a, d
        inc a ; 1,3,5
        ld bc, 0xFFFD 
        out (c), a ; coarse tune for channel (1,3,5)
        ld d, a
        ld a, (hl)
        and 0xf
        ld bc, 0xBFFD
        out (c), a
        ld a, e ; 1,2,3
        add a, 7 ; 8,9,10
        ld bc, 0xFFFD 
        out (c), a ; amplitude reg for channel (8,9,10)
        ld a, (hl)
        rra 
        rra 
        rra 
        rra 
        and 0xf
        ld bc, 0xBFFD
        out (c), a
        inc hl

        pop bc
        djnz ayloop

        ld a, (data+1) ; number of frames to skip between writes
        ld b, a
haltloop:        
        halt ; wait until next frame
        djnz haltloop

        ; check if the data has looped
        ld de, (datacount)
        dec de
        ld a, d
        or a, e
        jr nz, noloop
        ld de, (data+2) 
        ld hl, data+4
noloop:
        ld (datacount), de

        jp mainloop
        


data:
        INCBIN "../aydata.dat"
datacount:
        dw 0
    	
    	SAVETAP "test.tap", $8000   ; Save the assembled program as a tap file