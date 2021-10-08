//
//  DialPadRow.swift
//  Dial Pad
//
//  Created by dato che on 11/12/20.
//

import  UIKit


class DialPadRow: BaseReusableView{
    @IBOutlet var dialPadButton1: DialPadButton!
    @IBOutlet var dialPadButton2: DialPadButton!
    @IBOutlet var dialPadButton3: DialPadButton!
    
    weak var delegate: DialPadButtonDelegate?
    
    override func setup() {
        guard let delegate = delegate else { fatalError("delegate not set!") }
        dialPadButton1?.delegate = delegate
        dialPadButton2?.delegate = delegate
        dialPadButton3?.delegate = delegate
    }
    
    
}
