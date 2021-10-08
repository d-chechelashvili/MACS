//
//  ViewController.swift
//  Weather App
//
//  Created by dato che on 1/30/21.
//

import UIKit
import UPCarouselFlowLayout
import CoreLocation

class TodayViewController: UIViewController {

    @IBOutlet var carouselView: UICollectionView!
    @IBOutlet var pageControl: UIPageControl!
    
    @IBOutlet var errorPage: UIView!
    @IBOutlet var reloadButton: UIButton!
    
    @IBOutlet var activityIndicatorPage: UIView!
    @IBOutlet var activityIndicator: UIActivityIndicatorView!
    
    private var allCities = [String]()
    private var weathers = [ServiceResponse]()
    private var colors = [UIColor]()
    private let service = Service(for: .today)
    
    private let locationManager = CLLocationManager()
    
    
    override func viewDidLayoutSubviews() {
        super.viewDidLayoutSubviews()
        reloadButton.layer.cornerRadius = reloadButton.frame.width / 10
        setCarouselLayout()
    }
    
    private func setCarouselLayout() {
        let layout = UPCarouselFlowLayout()
        layout.scrollDirection = .horizontal
        layout.itemSize = CGSize(
            width: carouselView.frame.size.width * 0.74,
            height: carouselView.frame.size.height)
        //layout.spacingMode = .fixed(spacing: 15)
        carouselView.collectionViewLayout = layout
    }
    
    override func viewDidLoad() {
        super.viewDidLoad()
        configureCollectionView()
        makeBarsTransparent()
        loadCities()
        setupLocationManager()
    }
    
    private func configureCollectionView() {
        carouselView.delegate = self
        carouselView.dataSource = self
        carouselView.register(
            UINib(nibName: "TodayWeatherCell", bundle: nil),
            forCellWithReuseIdentifier: "weatherCell"
        )
        
        carouselView.addGestureRecognizer(
            UILongPressGestureRecognizer(
                target: self,
                action: #selector(handleLongPress(gesture:))
            )
        )
    }
    
    private func makeBarsTransparent() {
        self.navigationController?.navigationBar.setBackgroundImage(UIImage(), for: .default)
        self.navigationController?.navigationBar.shadowImage = UIImage()
        self.navigationController?.navigationBar.layoutIfNeeded()

        self.tabBarController?.tabBar.layer.borderWidth = 0
        self.tabBarController?.tabBar.clipsToBounds = true
        self.tabBarController?.tabBar.backgroundImage = UIImage()
        self.tabBarController?.tabBar.shadowImage = UIImage()
        self.tabBarController?.tabBar.layoutIfNeeded()
        self.tabBarController?.tabBar.unselectedItemTintColor = UIColor.white
    }
    
    private func setupLocationManager(){
        if CLLocationManager.locationServicesEnabled(){
            locationManager.delegate = self
            locationManager.desiredAccuracy = kCLLocationAccuracyNearestTenMeters
            locationManager.requestWhenInUseAuthorization()
            locationManager.startUpdatingLocation()
        }
    }
    
    private var group = DispatchGroup()
    private var errorHappened: Bool = false
    
    private func loadCities(){
        errorPage.isHidden = true
        activityIndicatorPage.isHidden = false
        allCities = UserDefaults.standard.array(forKey: "city.names") as? [String] ?? []
        for city in allCities{
            colors.append(.random)
            loadCity(city)
        }
        notifyDataLoad()
    }
    
    private func loadCity(_ city: String){
        group.enter()
        service.loadWeather(for: city){ [weak self] result in
            guard let self = self else {return}
            DispatchQueue.main.async {
                switch result {
                case .success(let weather):
                    self.weathers.append(weather)
                case .failure(_):
                    self.errorHappened = true
                }
                self.group.leave()
            }
        }
    }
    
    private func notifyDataLoad(){
        group.notify(queue: .main){
            self.activityIndicatorPage.isHidden = true
            if self.errorHappened{
                self.errorPage.isHidden = false
            }else{
                self.pageControl.numberOfPages = self.allCities.count
                self.carouselView.reloadData()
            }
            self.errorHappened = false
        }
    }
    
    
    @IBAction func handleAddCity(){
        let mainStoryBoard = UIStoryboard(name: "Main", bundle: nil)
        let vc = mainStoryBoard.instantiateViewController(withIdentifier: "AddCityViewController")
        vc.modalPresentationStyle = .overFullScreen
        self.present(vc, animated: true, completion: nil)
        
        if let addCityVC = vc as? AddCityViewController{
            addCityVC.delegate = self
        }
        
    }
    
    @IBAction func handleRefresh(){
        loadCities()
    }
    
    @objc func handleLongPress(gesture: UILongPressGestureRecognizer) {
        let location = gesture.location(in: carouselView)
        if let indexPath = carouselView.indexPathForItem(at: location) {
            if gesture.state == .began{
                deleteCity(indexPath)
            }
        }
    }
    
    private func deleteCity(_ indexPath: IndexPath){
        let cityName = allCities[indexPath.row]
        let messageText = "City " + cityName + " will be deleted"
        let alert = UIAlertController(title: "Delete City", message: messageText, preferredStyle: .alert)
        
        alert.addAction(
            UIAlertAction(
                title: "Cancel",
                style: .cancel,
                handler: nil
            )
        )
        
        alert.addAction(
            UIAlertAction(
                title: "Delete",
                style: .destructive,
                handler: { [unowned self] _ in
                    if let index = allCities.firstIndex(of: cityName){
                        allCities.remove(at: index)
                    }
                    UserDefaults.standard.set(allCities, forKey: "city.names")
                    if let index = weathers.firstIndex(where: {$0.name == cityName}){
                        weathers.remove(at: index)
                    }
                    pageControl.numberOfPages = allCities.count
                    carouselView.reloadData()
                }
            )
        )
        present(alert, animated: true, completion: nil)
    }
}
    
extension TodayViewController: UICollectionViewDelegate {
    func scrollViewDidScroll(_ scrollView: UIScrollView) {
        let pageWidth = scrollView.frame.size.width * 0.74
        let page = floor((scrollView.contentOffset.x - pageWidth/2)/pageWidth) + 1
        pageControl.currentPage = Int(page)
    }
}

extension TodayViewController: UICollectionViewDataSource {
    func numberOfSections(in collectionView: UICollectionView) -> Int {
        return 1
    }

    func collectionView(_ collectionView: UICollectionView, numberOfItemsInSection section: Int) -> Int{
        return allCities.count
    }
    
    func collectionView(_ collectionView: UICollectionView, cellForItemAt indexPath: IndexPath) -> UICollectionViewCell {
        let cell = initTodayWeatherCell(collectionView, indexPath)
        return cell
    }
    
    func initTodayWeatherCell(_ collectionView: UICollectionView,_ indexPath: IndexPath) -> UICollectionViewCell{
        let cell = collectionView.dequeueReusableCell(withReuseIdentifier: "weatherCell", for: indexPath)
        if let todayWeatherCell = cell as? TodayWeatherCell{
            todayWeatherCell.mainView.backgroundColor = colors[indexPath.row]
            let cityName = allCities[indexPath.row]
            if let index = weathers.firstIndex(where: {$0.name == cityName}){
                let weather = weathers[index]
                todayWeatherCell.cityLabel.text = weather.name + ", " + weather.sys.country
                todayWeatherCell.cloudnessLabel.text = weather.clouds.all.description + " %"
                todayWeatherCell.humidityLabel.text = weather.main.humidity.description + " mm"
                todayWeatherCell.weatherLabel.text = Int(round(weather.main.temp - 273.15)).description + "°C | " + weather.weather[0].main
                todayWeatherCell.windDirectionLabel.text = self.getWindDirection(weather.wind.deg)
                todayWeatherCell.windSpeedLabel.text = String(format: "%.2f",  (weather.wind.speed * 3.6)) + " km/h"
                todayWeatherCell.weatherImage.image = self.getWeatherImage(weather.weather[0].icon)
            }
        }
        return cell
    }
    
    func getWindDirection(_ degrees: Int) -> String{
        if degrees > 303 || degrees <= 56 {
           return "N"
        }
        else if degrees > 56 && degrees <= 123 {
            return "E"
        }
        else if degrees > 123 && degrees <= 236 {
            return "S"
        }
        return "W"
    }
    
    func getWeatherImage(_ icon: String) -> UIImage{
        let imageURl = "https://openweathermap.org/img/w/" + icon + ".png"
        guard let url = URL(string: imageURl) else { return UIImage() }
        if let data = try? Data(contentsOf: url){
            if let image = UIImage(data: data){
                return image
            }
        }
        return UIImage()
    }
    
}

extension TodayViewController: AddCityViewControllerDelegate{
    func addCityButtonTapped(_ sender: AddCityViewController, _ serviceResponse: ServiceResponse) {
        var cityNames: [String] = UserDefaults.standard.array(forKey: "city.names") as? [String] ?? []
        let cityName = serviceResponse.name
        if !cityNames.contains(cityName){
            cityNames.append(cityName)
            colors.append(.random)
            weathers.append(serviceResponse)
            UserDefaults.standard.set(cityNames, forKey: "city.names")
            allCities = cityNames
            pageControl.numberOfPages = allCities.count
            carouselView.reloadData()
        }
    }
}

extension TodayViewController: CLLocationManagerDelegate{
    func locationManager(_ manager: CLLocationManager, didUpdateLocations locations: [CLLocation]){
        guard let latitude = locationManager.location?.coordinate.latitude else { return  }
        guard let longtitude = locationManager.location?.coordinate.longitude else{ return }
        service.loadWeather(for: latitude, longtitude: longtitude){ [weak self] result in
            guard let self = self else {return}
            DispatchQueue.main.async {
                switch result {
                case .success(let weather):
                    var cityNames: [String] = UserDefaults.standard.array(forKey: "city.names") as? [String] ?? []
                    let cityName = weather.name
                    if !cityNames.contains(cityName){
                        cityNames.append(cityName)
                        self.colors.append(.random)
                        self.weathers.append(weather)
                        UserDefaults.standard.set(cityNames, forKey: "city.names")
                        self.allCities = cityNames
                        self.pageControl.numberOfPages = self.allCities.count
                            if let index = self.allCities.firstIndex(of: cityName){
                                DispatchQueue.main.async {
                                self.carouselView.scrollToItem(at: IndexPath(row: index, section: 0), at: .centeredHorizontally, animated: false)
                                }
                                self.pageControl.currentPage = index
                            }
                        self.carouselView.reloadData()
                    }else{
                        if let index = self.allCities.firstIndex(of: cityName){
                            DispatchQueue.main.async {
                            self.carouselView.scrollToItem(at: IndexPath(row: index, section: 0), at: .centeredHorizontally, animated: false)
                            }
                            self.pageControl.currentPage = index
                        }
                    }
                case .failure(_):
                    return
                }
            }
        }
    }
}

extension TodayViewController{
    func getSelectedCity() -> String?{
        if allCities.count == 0 {return nil}
        return allCities[pageControl.currentPage]
    }
}



extension UIColor {
    

    class var random: UIColor {
        return UIColor(
            red: .random(in: 0 ... 1),
            green: .random(in: 0 ... 1),
            blue: .random(in: 0 ... 1),
            alpha: 1
        )
    }
}
