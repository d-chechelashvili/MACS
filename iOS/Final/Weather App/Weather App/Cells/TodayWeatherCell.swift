//
//  TodayWeatherCell.swift
//  Weather App
//
//  Created by dato che on 2/3/21.
//

import UIKit

class TodayWeatherCell: UICollectionViewCell {
    
    @IBOutlet var mainView: UIView!
    @IBOutlet var weatherImage: UIImageView!
    @IBOutlet var cityLabel: UILabel!
    @IBOutlet var weatherLabel: UILabel!
    @IBOutlet var cloudnessLabel: UILabel!
    @IBOutlet var humidityLabel: UILabel!
    @IBOutlet var windSpeedLabel: UILabel!
    @IBOutlet var windDirectionLabel: UILabel!
    
    override func awakeFromNib() {
        super.awakeFromNib()
        mainView.layer.cornerRadius = mainView.frame.height / 32
    }

}
