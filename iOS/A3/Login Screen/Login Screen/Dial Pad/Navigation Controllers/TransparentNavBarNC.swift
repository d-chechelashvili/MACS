//
//  TransparentNavBarNC.swift
//  Dial Pad
//
//  Created by dato che on 11/13/20.
//

import UIKit

class TransparentNavBarNC: UINavigationController {
    
    override func viewDidLoad() {
        super.viewDidLoad()
    
        let appearece = UINavigationBarAppearance()
        appearece.configureWithTransparentBackground()
        
        navigationBar.compactAppearance = appearece
        navigationBar.standardAppearance = appearece
        navigationBar.scrollEdgeAppearance = appearece
    }
    
}
